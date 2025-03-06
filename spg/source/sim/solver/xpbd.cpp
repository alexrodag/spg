#include <spg/sim/solver/xpbd.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/sim/rigidBodyGroup.h>
#include <spg/utils/timer.h>
#include <spg/utils/graphColoring.h>
#include <spg/utils/functionalUtilities.h>

#include <iostream>

namespace spg::solver
{
XPBD::XPBD(const bool useParallelGaussSeidel)
{
    m_useParallelGaussSeidel = useParallelGaussSeidel;
}
void XPBD::step()
{
    if (m_useParallelGaussSeidel && m_parallelGroupsDirty) {
        computeParallelStencilGroups();
        m_parallelGroupsDirty = false;
    }
    const Real dt = m_dtStep / m_nsubsteps;
    const int nObjects = static_cast<int>(std::get<std::vector<SimObject>>(m_objects).size());
    m_simObjectsOldPos.resize(nObjects);
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "XPBD step\n";
    }
    Timer timer;
    timer.start();
    Timer detailTimer;
    Real inertialTime = 0;
    Real projectionTime = 0;
    Real updateTime = 0;
    for (int s = 0; s < m_nsubsteps; ++s) {
        // Compute predicted inertial pos
        detailTimer.start();
        for (int objId = 0; objId < nObjects; ++objId) {
            auto &object = std::get<std::vector<SimObject>>(m_objects)[objId];
            auto &positions = object.positions();
            m_simObjectsOldPos[objId] = positions;
            const auto &velocities = object.velocities();
            const auto &invMass = object.invMasses();
            const int nParticles = static_cast<int>(object.nElements());
            const Vector3 dtdtg = dt * dt * m_gravity;
            for (int i = 0; i < nParticles; ++i) {
                if (invMass[i] != 0) {
                    positions[i] = positions[i] + dt * velocities[i] + dtdtg;
                }
            }
        }
        detailTimer.stop();
        inertialTime += detailTimer.getMilliseconds();

        // Project constraints positions
        detailTimer.start();
        for (auto &obj : std::get<std::vector<SimObject>>(m_objects)) {
            const auto &energies = obj.energies();
            for (const auto &energy : energies) {
                if (auto stencilGroupsIt = m_stencilGroupsPerEnergy.find(static_cast<void *>(energy.get()));
                    stencilGroupsIt != m_stencilGroupsPerEnergy.end()) {
                    // If we have computed stencil groups for the energy, use Parallel Gauss Seidel path, inspired by
                    // ref "Vivace: a Practical Gauss-Seidel Method for Stable Soft Body Dynamics"
                    for (const auto &group : stencilGroupsIt->second) {
                        const int groupSize = static_cast<int>(group.size());
#pragma omp parallel for
                        for (int i = 0; i < groupSize; ++i) {
                            energy->projectPosition(group[i], obj, dt);
                        }
                    }
                } else {
                    // Serial Gauss Seidel path otherwise
                    const auto nstencils = energy->nStencils();
                    for (int i = 0; i < nstencils; ++i) {
                        energy->projectPosition(i, obj, dt);
                    }
                }
            }
        }
        detailTimer.stop();
        projectionTime += detailTimer.getMilliseconds();

        // Update velocities
        detailTimer.start();
        for (int objId = 0; objId < nObjects; ++objId) {
            auto &object = std::get<std::vector<SimObject>>(m_objects)[objId];
            const auto &positions = object.positions();
            const auto &oldPositions = m_simObjectsOldPos[objId];
            auto &velocities = object.velocities();
            const int nParticles = static_cast<int>(object.nElements());
            for (int i = 0; i < nParticles; ++i) {
                velocities[i] = (positions[i] - oldPositions[i]) / dt;
            }
            // basic damping
            // TODO: Change for some unified damping across solvers
            /* for (int i = 0; i < nParticles; ++i) {
                velocities[i] = velocities[i] * 0.999;
            } */
        }
        detailTimer.stop();
        updateTime += detailTimer.getMilliseconds();
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Inertial position time: " << inertialTime << "ms\n";
        std::cout << "  Projection time: " << projectionTime << "ms\n";
        std::cout << "  Update time: " << updateTime << "ms\n";
        std::cout << "  Total step time: " << timer.getMilliseconds() << "ms\n\n";
    }
}

void XPBD::computeParallelStencilGroups()
{
    // TODO: Here we are coloring all energies for the sake of simplicity. Performance-wise it is generally worth
    // leaving energies with low number of stencils to compute in a serial manner, as the overhead of thread managing is
    // very likely to outweight the gain due to parallel computation.
    Timer timer;
    timer.start();
    m_stencilGroupsPerEnergy.clear();
    apply_each(
        [this](const auto &objs) {
            for (const auto &obj : objs) {
                const auto &energies = obj.energies();
                for (const auto &energy : energies) {
                    const auto stencilColors = coloring::colorStencils({energy->flatStencils(), energy->stencilSize()});
                    auto &stencilGroups = m_stencilGroupsPerEnergy[energy.get()];
                    for (int i = 0; i < stencilColors.size(); ++i) {
                        stencilGroups.resize(std::max(stencilGroups.size(), static_cast<size_t>(stencilColors[i] + 1)));
                        stencilGroups[stencilColors[i]].push_back(i);
                    }
                }
            }
        },
        m_objects);
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "XPBD stencil coloring time: " << timer.getMilliseconds() << "ms\n";
    }
}
}  // namespace spg::solver