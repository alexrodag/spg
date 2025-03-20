#include <spg/sim/solver/vbd.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject/particleGroup.h>
#include <spg/sim/simObject/rigidBodyGroup.h>
#include <spg/utils/timer.h>
#include <spg/utils/graphColoring.h>
#include <spg/utils/functionalUtilities.h>

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
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "VDB step\n";
    }
    Timer timer;
    timer.start();
    Timer detailTimer;
    Real inertialTime = 0;
    Real descentTime = 0;
    Real updateTime = 0;

    // Compute total DOFs
    int accumulatedNDOF = 0;
    apply_each(
        [&accumulatedNDOF](const auto &objs) {
            for (const auto &obj : objs) {
                accumulatedNDOF += obj.nDOF();
            }
        },
        m_objects);
    const int totalNDOF{accumulatedNDOF};
    m_xOld.resize(totalNDOF);
    m_xInertial.resize(totalNDOF);

    for (int s = 0; s < m_nsubsteps; ++s) {
        detailTimer.start();
        // Store previous positions
        apply_each(
            [this](const auto &objs) {
                int accumulatedNDOF = 0;
                for (const auto &obj : objs) {
                    obj.getPositions(m_xOld, accumulatedNDOF);
                    accumulatedNDOF += obj.nDOF();
                }
            },
            m_objects);
        // Apply explicit forces to the velocities and compute predicted inertial pos and initial guess
        apply_each(
            [this, dt](auto &objs) {
                int accumulatedNDOF = 0;
                const Vector3 dtg = dt * m_gravity;
                auto l_skew = [](const spg::Vector3 &v) {
                    spg::Matrix3 vSkew;
                    vSkew << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
                    return vSkew;
                };
                for (auto &obj : objs) {
                    const int nPrimitives = static_cast<int>(obj.size());
                    for (int i = 0; i < nPrimitives; ++i) {
                        // TODO: This can be made more efficient by merging operations and bypassing the velocity, but I
                        // leave it as it is for the sake of clarity
                        obj.velocities()[i] += dtg;
                        if constexpr (std::is_same_v<std::decay_t<decltype(obj)>, RigidBodyGroup>) {
                            obj.omegas()[i] += dt * obj.invInertias()[i] *
                                               (-l_skew(obj.omegas()[i]) * obj.inertias()[i] * obj.omegas()[i]);
                        }
                    }
                    obj.integrateVelocities(dt);
                    obj.getPositions(m_xInertial, accumulatedNDOF);
                    accumulatedNDOF += obj.nDOF();
                }
                if (m_initialGuessType == InitialGuessType::InertialWithAcceleration) {
                    // Guess with constant velocity and explicit accelerations, already stored in objects state, no need
                    // to do anything
                }
            },
            m_objects);
        detailTimer.stop();
        inertialTime += detailTimer.getMilliseconds();

        // Iterative per node coordinate block descent
        detailTimer.start();
        const int iterations = m_iterations;
        for (int iter = 0; iter < iterations; ++iter) {
            int objId = 0;
            int accumulatedNDOF = 0;
            apply_each(
                [this, &objId, &accumulatedNDOF, invdt, invdtSquared, epsilon, iter](auto &objs) {
                    for (auto &obj : objs) {
                        auto &elementsPerVertex = m_simObjectsElementsPerVertex[objId];
                        const int nVertices = obj.size();
                        // Vertex descent function
                        auto l_vertexDescent =
                            [this, &obj, &elementsPerVertex, accumulatedNDOF, invdt, invdtSquared, epsilon, iter](
                                const int vIdx) {
                                const auto mOverdtSquared = obj.masses()[vIdx] * invdtSquared;
                                Vector<obj.s_nDOFs> f;
                                Matrix<obj.s_nDOFs, obj.s_nDOFs> H;
                                f.setZero();
                                H.setZero();
                                if (iter != 0 || m_initialGuessType == InitialGuessType::InertialWithAcceleration) {
                                    // The first iter with inertial guess produces no forces here, so we can skip it
                                    f.segment<3>(0) = -mOverdtSquared *
                                                      (obj.positions()[vIdx] -
                                                       m_xInertial.segment<3>(vIdx * obj.s_nDOFs + accumulatedNDOF));
                                    if constexpr (std::is_same_v<std::decay_t<decltype(obj)>, RigidBodyGroup>) {
                                        f.segment<3>(3) =
                                            -invdt * (obj.inertias()[vIdx] *
                                                      (obj.computeIntegratedOmega(
                                                           m_xOld.segment<3>(vIdx * obj.s_nDOFs + 3 + accumulatedNDOF),
                                                           invdt,
                                                           vIdx) -
                                                       obj.omegas()[vIdx]));
                                    }
                                }
                                H.block<3, 3>(0, 0) = mOverdtSquared * Matrix3::Identity();
                                if constexpr (std::is_same_v<std::decay_t<decltype(obj)>, RigidBodyGroup>) {
                                    H.block<3, 3>(3, 3) = invdtSquared * obj.inertias()[vIdx];
                                }
                                auto &vertexElementEntriesPerEnergy = elementsPerVertex[vIdx];
                                for (int e = 0; e < obj.energies().size(); ++e) {
                                    auto energy = obj.energies()[e].get();
                                    const auto &vertexElementEntries = vertexElementEntriesPerEnergy[e];
                                    // Note: This inner loop would be implemented in parallel in the GPU version,
                                    // but for multicore CPU it does not improve performance, as the outer parallel
                                    // loop already feeds all the cores in current CPUs
                                    for (const auto [stencilIdx, vertexIdx] : vertexElementEntries) {
                                        energy->accumulateVertexForce(stencilIdx, vertexIdx, obj, f);
                                        energy->accumulateVertexHessian(stencilIdx, vertexIdx, obj, H);
                                    }
                                }
                                if (H.determinant() > epsilon) {
                                    // Note: If required, it is possible to do a PSD projection to the Hessian here,
                                    // before inverting it
                                    const auto deltax = (H.inverse() * f).eval();
                                    // Note: Line search could be added here as well, but they mention in the paper
                                    // that it is typically not needed
                                    obj.updateElementPositionFromDx(deltax, vIdx);
                                }
                            };
                        if (!m_simObjectsVertexGroups.empty()) {
                            // If vertex groups have been computed, run vertex descent in a Parallel Gauss Seidel
                            // fashion
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
                        accumulatedNDOF += obj.nDOF();
                    }
                },
                m_objects);
        }
        detailTimer.stop();
        descentTime += detailTimer.getMilliseconds();

        // Update velocities
        detailTimer.start();
        accumulatedNDOF = 0;
        apply_each(
            [this, &accumulatedNDOF, invdt](auto &objs) {
                for (auto &obj : objs) {
                    obj.computeIntegratedVelocities(m_xOld, accumulatedNDOF, invdt);
                    accumulatedNDOF += obj.nDOF();
                }
            },
            m_objects);
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
}

void VBD::computeParallelVertexGroups()
{
    Timer timer;
    timer.start();
    m_simObjectsVertexGroups.clear();
    m_simObjectsVertexGroups.resize(numSimObjects());
    int objId = 0;
    apply_each(
        [this, &objId](const auto &objs) {
            for (const auto &obj : objs) {
                auto &vertexGroups = m_simObjectsVertexGroups[objId];
                objId++;
                std::vector<coloring::FlatStencils> flatStencilsSet;
                for (const auto &energy : obj.energies()) {
                    flatStencilsSet.push_back({energy->flatStencils(), energy->stencilSize()});
                }
                auto vertexColors = coloring::colorVertices(obj.size(), flatStencilsSet);
                for (int i = 0; i < vertexColors.size(); ++i) {
                    vertexGroups.resize(std::max(vertexGroups.size(), static_cast<size_t>(vertexColors[i] + 1)));
                    vertexGroups[vertexColors[i]].push_back(i);
                }
            }
        },
        m_objects);
    timer.stop();
    std::cout << "VBD vertex coloring time: " << timer.getMilliseconds() << "ms\n";
}

void VBD::computeStencilInfoPerVertex()
{
    Timer timer;
    timer.start();
    m_simObjectsElementsPerVertex.resize(numSimObjects());
    int objId = 0;
    apply_each(
        [this, &objId](const auto &objs) {
            for (const auto &obj : objs) {
                const int nVertices = obj.size();
                const auto &energies = obj.energies();
                const int nEnergies = static_cast<int>(energies.size());
                auto &elementsPerVertex = m_simObjectsElementsPerVertex[objId];
                objId++;
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
        },
        m_objects);
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "VBD stencil incidence per vertex computation time: " << timer.getMilliseconds() << "ms\n";
    }
}
}  // namespace spg::solver