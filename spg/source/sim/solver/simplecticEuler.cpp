#include <spg/sim/solver/simplecticEuler.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/utils/timer.h>

#include <iostream>

namespace spg::solver
{
void SimplecticEuler::step()
{
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "SimplecticEuler step\n";
    }
    const Real dt = m_dtStep / m_nsubsteps;
    const int nObjects = static_cast<int>(m_simObjects.size());
    Timer timer;
    timer.start();
    for (int s = 0; s < m_nsubsteps; ++s) {
        // Update velocities
        for (int objId = 0; objId < nObjects; ++objId) {
            auto &object = m_simObjects[objId];
            const auto &energies = object.energies();
            for (auto &energy : energies) {
                const auto nstencils = energy->nStencils();
#pragma omp parallel for
                for (int i = 0; i < nstencils; ++i) {
                    energy->updateExplicitVelocities(i, object, dt, true);
                }
            }
            const auto &invMass = object.invMasses();
            auto &velocities = object.velocities();
            const int nParticles = static_cast<int>(object.nParticles());
            const Vector3 dtg = dt * m_gravity;
            for (int i = 0; i < nParticles; ++i) {
                if (invMass[i] != 0) {
                    velocities[i] += dtg;
                }
            }
        }

        // Update positions
        for (int objId = 0; objId < nObjects; ++objId) {
            auto &object = m_simObjects[objId];
            auto &positions = object.positions();
            const auto &velocities = object.velocities();
            const auto &invMass = object.invMasses();
            const int nParticles = static_cast<int>(object.nParticles());
            for (int i = 0; i < nParticles; ++i) {
                if (invMass[i] != 0) {
                    positions[i] = positions[i] + dt * velocities[i];
                }
            }
        }
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Total step time: " << timer.getMilliseconds() << "ms\n\n";
    }
}
}  // namespace spg::solver