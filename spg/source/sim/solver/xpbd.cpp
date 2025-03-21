#include <spg/sim/solver/xpbd.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject/particleGroup.h>
#include <spg/sim/simObject/rigidBodyGroup.h>
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
    const Real invdt = 1 / dt;
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
        detailTimer.start();
        accumulatedNDOF = 0;
        // Store previous positions
        apply_each(
            [this, &accumulatedNDOF](const auto &objs) {
                for (const auto &obj : objs) {
                    obj.getPositions(m_xOld, accumulatedNDOF);
                    accumulatedNDOF += obj.nDOF();
                }
            },
            m_objects);
        // Apply explicit forces to the velocities and compute predicted inertial pos
        apply_each(
            [this, dt](auto &objs) {
                auto l_skew = [](const spg::Vector3 &v) {
                    spg::Matrix3 vSkew;
                    vSkew << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
                    return vSkew;
                };
                const Vector3 dtg = dt * m_gravity;
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
                }
            },
            m_objects);
        detailTimer.stop();
        inertialTime += detailTimer.getMilliseconds();

        // Project constraints positions
        detailTimer.start();
        apply_each(
            [this, dt](auto &objs) {
                for (auto &obj : objs) {
                    const auto &energies = obj.energies();
                    for (const auto &energy : energies) {
                        if (auto stencilGroupsIt = m_stencilGroupsPerEnergy.find(static_cast<void *>(energy.get()));
                            stencilGroupsIt != m_stencilGroupsPerEnergy.end()) {
                            // If we have computed stencil groups for the energy, use Parallel Gauss Seidel path,
                            // inspired by ref "Vivace: a Practical Gauss-Seidel Method for Stable Soft Body Dynamics"
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
            },
            m_objects);
        detailTimer.stop();
        projectionTime += detailTimer.getMilliseconds();

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