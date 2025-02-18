#include <spg/sim/solver/vbd.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/sim/rigidBodyGroup.h>
#include <spg/utils/timer.h>
#include <spg/utils/graphColoring.h>

#include <iostream>

namespace spg::solver
{
void VBD::step()
{
    if (m_precomputationUpdateRequired) {
        computeParallelVertexGroups();
        computeStencilInfoPerVertex();
        m_precomputationUpdateRequired = false;
    }
    const Real dt = m_dtStep / m_nsubsteps;
    const Real invdt = 1. / dt;
    const Real invdtSquared = invdt * invdt;
    constexpr Real epsilon = 1e-10;
    const int nObjects = static_cast<int>(std::get<std::vector<SimObject>>(m_objects).size());
    m_simObjectsOldPos.resize(nObjects);
    m_simObjectsInertialPositions.resize(nObjects);
    m_simObjPrevStepVelocities.resize(nObjects);
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "VDB step\n";
    }
    Timer timer;
    timer.start();
    Timer detailTimer;
    Real inertialTime = 0;
    Real descentTime = 0;
    Real updateTime = 0;

    for (int s = 0; s < m_nsubsteps; ++s) {
        // Compute predicted inertial pos and initial guess per vertex
        detailTimer.start();
        for (int objId = 0; objId < nObjects; ++objId) {
            auto &object = std::get<std::vector<SimObject>>(m_objects)[objId];
            auto &positions = object.positions();
            const auto &velocities = object.velocities();
            const int nVertices = object.nElements();
            const Real dtdt = dt * dt;
            const Vector3 dtdtg = dtdt * m_gravity;
            m_simObjectsOldPos[objId] = positions;
            auto &inertialPositions = m_simObjectsInertialPositions[objId];
            inertialPositions.resize(nVertices);
            // Set prevStep velocities with current ones if first time or update required
            auto &prevStepVelocities = m_simObjPrevStepVelocities[objId];
            if (prevStepVelocities.size() != velocities.size()) {
                prevStepVelocities = velocities;
            }
            for (int i = 0; i < nVertices; ++i) {
                inertialPositions[i] = positions[i] + dt * velocities[i] + dtdtg;
            }
            if (m_initialGuessType == InitialGuessType::Inertial) {
                // Guess with constant velocity
                positions = inertialPositions;
                std::for_each(positions.begin(), positions.end(), [dtdtg](Vector3 &pos) { pos -= dtdtg; });
            } else if (m_initialGuessType == InitialGuessType::InertialWithAcceleration) {
                // Guess with constant velocity and explicit accelerations
                positions = inertialPositions;
            } else if (m_initialGuessType == InitialGuessType::Adaptive) {
                // Guess with constant velocity and adaptive accelerations based on the paper
                // Note: Doesnt seem to work that well, introduces oscillation. May be worth investigating why
                const Real aExtNorm = m_gravity.norm();
                const Vector3 aExtDirection = m_gravity / aExtNorm;
                std::vector<Vector3> adaptiveAcceleration(nVertices);
                // TODO: This loop has potential for parallelism, profile.
                for (int i = 0; i < nVertices; ++i) {
                    adaptiveAcceleration[i] = (velocities[i] - prevStepVelocities[i]) / dt;
                    const Real aExtProjection = adaptiveAcceleration[i].dot(aExtDirection);
                    const Real factor =
                        std::clamp(aExtProjection / aExtNorm, static_cast<Real>(0), static_cast<Real>(1));
                    adaptiveAcceleration[i] = m_gravity * factor;
                    positions[i] = positions[i] + dt * velocities[i] + adaptiveAcceleration[i] * dtdt;
                }
            }
            // Store prev velocities for next step
            m_simObjPrevStepVelocities[objId] = velocities;
        }
        detailTimer.stop();
        inertialTime += detailTimer.getMilliseconds();

        // Iterative per node coordinate block descent
        detailTimer.start();
        const int iterations = m_iterations;
        for (int iter = 0; iter < iterations; ++iter) {
            for (int objId = 0; objId < std::get<std::vector<SimObject>>(m_objects).size(); ++objId) {
                auto &obj = std::get<std::vector<SimObject>>(m_objects)[objId];
                const auto &inertialPositions = m_simObjectsInertialPositions[objId];
                auto &elementsPerVertex = m_simObjectsElementsPerVertex[objId];
                const int nVertices = obj.nElements();
                // Vertex descent function
                auto l_vertexDescent =
                    [&obj, &inertialPositions, &elementsPerVertex, invdtSquared, epsilon](const int vIdx) {
                        const auto mOverdtSquared = obj.masses()[vIdx] * invdtSquared;
                        Vector3 f = -mOverdtSquared * (obj.positions()[vIdx] - inertialPositions[vIdx]);
                        Matrix3 H = mOverdtSquared * Matrix3::Identity();
                        auto &vertexElementEntriesPerEnergy = elementsPerVertex[vIdx];
                        for (int e = 0; e < obj.energies().size(); ++e) {
                            auto energy = obj.energies()[e].get();
                            const auto &vertexElementEntries = vertexElementEntriesPerEnergy[e];
                            // Note: This inner loop would be implemented in parallel in the GPU version, but for
                            // multicore CPU it does not improve performance, as the outer parallel loop
                            // already feeds all the cores in current CPUs
                            for (const auto [stencilIdx, vertexIdx] : vertexElementEntries) {
                                energy->accumulateVertexForce(stencilIdx, vertexIdx, obj, f);
                                energy->accumulateVertexHessian(stencilIdx, vertexIdx, obj, H);
                            }
                        }
                        if (H.determinant() > epsilon) {
                            // Note: If required, it is possible to do a PSD projection to the Hessian here, before
                            // inverting it
                            const Vector3 deltax = H.inverse() * f;
                            // Note: Line search could be added here as well, but they mention in the paper that it is
                            // typically not needed
                            obj.positions()[vIdx] += deltax;
                        }
                    };
                if (!m_simObjectsVertexGroups.empty()) {
                    // If vertex groups have been computed, run vertex descent in a Parallel Gauss Seidel fashion
                    const auto &vertexGroups = m_simObjectsVertexGroups[objId];
                    for (const auto &vertexGroup : vertexGroups) {
                        const int verticesInGroup = static_cast<int>(vertexGroup.size());
#pragma omp parallel for
                        for (int v = 0; v < verticesInGroup; ++v) {
                            l_vertexDescent(vertexGroup[v]);
                        }
                    }
                } else {
                    // Run serial Gauss Seidel otherwise
                    for (int v = 0; v < nVertices; ++v) {
                        l_vertexDescent(v);
                    }
                }
            }
        }
        detailTimer.stop();
        descentTime += detailTimer.getMilliseconds();

        // Update velocities
        detailTimer.start();
        for (int objId = 0; objId < nObjects; ++objId) {
            auto &object = std::get<std::vector<SimObject>>(m_objects)[objId];
            const auto &positions = object.positions();
            const auto &oldPositions = m_simObjectsOldPos[objId];
            auto &velocities = object.velocities();
            const int nParticles = static_cast<int>(object.nElements());
            for (int i = 0; i < nParticles; ++i) {
                velocities[i] = (positions[i] - oldPositions[i]) * invdt;
            }
        }
        detailTimer.stop();
        updateTime += detailTimer.getMilliseconds();
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        const auto totalTime = timer.getMilliseconds();
        std::cout << "  Inertial position time: " << inertialTime << "ms\n";
        std::cout << "  Coordinate block descent time: " << descentTime << "ms\n";
        std::cout << "  Update time: " << updateTime << "ms\n";
        std::cout << "  Total step time: " << totalTime << "ms\n\n";
    }
}

void VBD::reset()
{
    BaseSolver::reset();
    m_precomputationUpdateRequired = true;
    m_simObjPrevStepVelocities.clear();
}

void VBD::requirePrecomputationUpdate()
{
    m_precomputationUpdateRequired = true;
}

void VBD::computeParallelVertexGroups()
{
    Timer timer;
    timer.start();
    m_simObjectsVertexGroups.clear();
    m_simObjectsVertexGroups.resize(std::get<std::vector<SimObject>>(m_objects).size());
    for (int objId = 0; objId < std::get<std::vector<SimObject>>(m_objects).size(); ++objId) {
        const auto &object = std::get<std::vector<SimObject>>(m_objects)[objId];
        auto &vertexGroups = m_simObjectsVertexGroups[objId];
        std::vector<coloring::FlatStencils> flatStencilsSet;
        for (const auto &energy : object.energies()) {
            flatStencilsSet.push_back({energy->flatStencils(), energy->stencilSize()});
        }
        auto vertexColors = coloring::colorVertices(object.nElements(), flatStencilsSet);
        for (int i = 0; i < vertexColors.size(); ++i) {
            vertexGroups.resize(std::max(vertexGroups.size(), static_cast<size_t>(vertexColors[i] + 1)));
            vertexGroups[vertexColors[i]].push_back(i);
        }
    }
    timer.stop();
    std::cout << "VBD vertex coloring time: " << timer.getMilliseconds() << "ms\n";
}

void VBD::computeStencilInfoPerVertex()
{
    Timer timer;
    timer.start();
    m_simObjectsElementsPerVertex.resize(std::get<std::vector<SimObject>>(m_objects).size());
    for (int objId = 0; objId < std::get<std::vector<SimObject>>(m_objects).size(); ++objId) {
        const auto &obj = std::get<std::vector<SimObject>>(m_objects)[objId];
        const int nVertices = obj.nElements();
        const auto &energies = obj.energies();
        const int nEnergies = static_cast<int>(energies.size());
        auto &elementsPerVertex = m_simObjectsElementsPerVertex[objId];
        elementsPerVertex.resize(nVertices);
        for (int v = 0; v < nVertices; ++v) {
            elementsPerVertex[v].clear();
            elementsPerVertex[v].resize(nEnergies);
        }
        for (int e = 0; e < energies.size(); ++e) {
            const auto &energy = energies[e];
            const int stencilSize = energy->stencilSize();
            const auto flatStencils = energy->flatStencils();
            const int nStencilEntries = static_cast<int>(flatStencils.size());
            for (int s = 0; s < nStencilEntries; ++s) {
                const int vertexIdx = flatStencils[s];
                const ElementEntry entry = {s / stencilSize, s % stencilSize};
                elementsPerVertex[vertexIdx][e].push_back(entry);
            }
        }
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "VBD stencil incidence per vertex computation time: " << timer.getMilliseconds() << "ms\n";
    }
}
}  // namespace spg::solver