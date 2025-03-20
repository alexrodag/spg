#include <spg/sim/solver/simplecticEuler.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject/particleGroup.h>
#include <spg/sim/simObject/rigidBodyGroup.h>
#include <spg/utils/timer.h>
#include <spg/utils/functionalUtilities.h>

#include <iostream>

namespace spg::solver
{
void SimplecticEuler::step()
{
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "SimplecticEuler step\n";
    }
    const Real dt = m_dtStep / m_nsubsteps;
    Timer timer;
    timer.start();
    for (int s = 0; s < m_nsubsteps; ++s) {
        // Update velocities
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
                    // Accumulate velocities due to external forces
                    for (int i = 0; i < nPrimitives; ++i) {
                        auto &velocities = obj.velocities();
                        velocities[i] += dtg;
                        if constexpr (std::is_same_v<std::decay_t<decltype(obj)>, RigidBodyGroup>) {
                            obj.omegas()[i] += dt * obj.invInertias()[i] *
                                               (-obj.omegas()[i].cross(obj.inertias()[i] * obj.omegas()[i]));
                        }
                    }
                    // Accumulate velocities due to external forces
                    const auto &energies = obj.energies();
                    for (auto &energy : energies) {
                        const auto nstencils = energy->nStencils();
#pragma omp parallel for
                        for (int i = 0; i < nstencils; ++i) {
                            energy->updateExplicitVelocities(i, obj, dt, true);
                        }
                    }
                }
            },
            m_objects);

        // Update positions
        apply_each(
            [dt](auto &objs) {
                for (auto &obj : objs) {
                    obj.integrateVelocities(dt);
                }
            },
            m_objects);
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Total step time: " << timer.getMilliseconds() << "ms\n\n";
    }
}
}  // namespace spg::solver