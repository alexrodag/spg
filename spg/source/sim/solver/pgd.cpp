#include <spg/sim/solver/pgd.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/utils/timer.h>

#include <iostream>

namespace spg::solver
{
void PGD::step()
{
    const Real dt = m_dtStep / m_nsubsteps;
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "PGD step\n";
    }
    Timer timer;
    timer.start();
    Timer detailTimer;
    Real forceTime = 0;
    Real hessianTime = 0;
    Real energyTime = 0;
    int accumulatedNDOF = 0;
    for (const auto &object : m_simObjects) {
        accumulatedNDOF += object.nDOF();
    }
    const int totalNDOF{accumulatedNDOF};
    VectorX x0(totalNDOF);
    VectorX v0(totalNDOF);
    VectorX xBackup(totalNDOF);
    VectorX vBackup(totalNDOF);
    VectorX vTilde(totalNDOF);
    VectorX v(totalNDOF);
    VectorX vCandidate(totalNDOF);
    VectorX f(totalNDOF);
    VectorX diagonalMass(totalNDOF);
    for (int s = 0; s < m_nsubsteps; ++s) {
        getSystemState(x0, v0);
        if (m_verbosity == Verbosity::Debug) {
            std::cout << "Starting step\n";
        }
        // Compute inertial velocities (including explicit accelerations) and update state
        vTilde = v0;
        const Vector3 dtg = dt * m_gravity;
        for (int i = 0; i < totalNDOF / 3; ++i) {
            vTilde.segment<3>(i * 3) += dtg;
        }
        // Set initial guess to be the inertial velocity with explicit accelerations.
        if (m_initialGuessType == InitialGuessType::Inertial) {
            // Guess with constant velocity
            v = v0;
        } else if (m_initialGuessType == InitialGuessType::InertialWithAcceleration) {
            // Guess with constant velocity and explicit accelerations.
            v = vTilde;
        } else if (m_initialGuessType == InitialGuessType::Adaptive) {
            // Guess with constant velocity and adaptive accelerations based on ref: "Vertex Block Descent".
            v = v0;
            VectorX aPrev(totalNDOF);
            if (m_vPrevStep.size() != v0.size()) {
                m_vPrevStep = v0;
            }
            aPrev = (v0 - m_vPrevStep) / dt;
            const Real aExtNorm = m_gravity.norm();
            const Vector3 aDirection = m_gravity / aExtNorm;
            for (int i = 0; i < totalNDOF / 3; ++i) {
                Real factor = std::clamp(
                    aPrev.segment<3>(i * 3).dot(aDirection) / aExtNorm, static_cast<Real>(0.0), static_cast<Real>(1.0));
                v.segment<3>(i * 3) += dtg * factor;
            }
        }

        VectorX x = x0 + dt * v;
        setObjectsState(x, v);
        getSystemDiagMass(diagonalMass);
        bool badConvergence = false;
        bool recomputePreconditioner = false;
        for (int i = 0; i < m_iterations; ++i) {
            if (m_verbosity == Verbosity::Debug) {
                std::cout << "Iteration " << i << "\n";
            }
            m_alpha = s_maxAlpha;
            xBackup = x;
            vBackup = v;
            if (i % m_preconditionerUpdateCadency == 0 || recomputePreconditioner) {
                if (m_verbosity == Verbosity::Debug) {
                    std::cout << "Recomputing preconditioner in iteration " << i << "\n";
                }
                // Compute diagonal of the system hessian
                m_preconditioner.resize(totalNDOF);
                detailTimer.start();
                getSystemDiagStiffness(m_preconditioner);
                detailTimer.stop();
                hessianTime += detailTimer.getMilliseconds();
                m_preconditioner *= -(dt * dt);
                m_preconditioner += diagonalMass;
                // Project diagonal hessian to absolute to ensure SPD
                m_preconditioner = m_preconditioner.cwiseAbs().eval();
                // Compute inverse of the hessian, to obtain the preconditioner
                m_preconditioner = m_preconditioner.cwiseInverse().eval();
            }
            detailTimer.start();
            getSystemForce(f);
            detailTimer.stop();
            forceTime += detailTimer.getMilliseconds();
            VectorX descentDirection = diagonalMass.cwiseProduct(v - vTilde) - dt * f;
            descentDirection = descentDirection.cwiseProduct(m_preconditioner).eval();

            // Line search
            // Note: This search takes a lot of time w.r.t. the rest of the computation. In "Descent Methods for Elastic
            // Body Simulation on the GPU" they suggest to do it spanning several interations (~8) to reduce its
            // performance impact, even if that means redoing more iterations when needed.
            detailTimer.start();
            const auto initialError = getSystemEnergy(v, vTilde, diagonalMass);
            detailTimer.stop();
            energyTime += detailTimer.getMilliseconds();
            if (m_verbosity == Verbosity::Debug) {
                std::cout << "  Pre iteration error: " << initialError << "\n";
            }
            auto error = initialError;
            constexpr Real marginFactor = 1.001;
            badConvergence = false;
            recomputePreconditioner = false;
            do {
                if (m_verbosity == Verbosity::Debug) {
                    std::cout << "  Doing line search with alpha " << m_alpha << "\n";
                }
                vCandidate = v - m_alpha * descentDirection;
                x = x0 + dt * vCandidate;
                setObjectsState(x, vCandidate);
                detailTimer.start();
                error = getSystemEnergy(vCandidate, vTilde, diagonalMass);
                detailTimer.stop();
                energyTime += detailTimer.getMilliseconds();
                if (m_verbosity == Verbosity::Debug) {
                    std::cout << "  New error: " << error << "\n";
                }
                m_alpha *= m_alphaSplittingFactor;
                if (m_alpha < s_maxAlpha * 0.5) {
                    recomputePreconditioner = true;
                }
                if (m_alpha < 0.01) {
                    badConvergence = true;
                }
            } while (error > (initialError * marginFactor) && !badConvergence);
            if (badConvergence) {
                m_alpha = s_maxAlpha;
                x = xBackup;
                v = vBackup;
                setObjectsState(x, v);
                if (m_verbosity == Verbosity::Debug) {
                    std::cout << "  Bad convergence, resetting to pre-iteration state\n";
                }
            } else {
                v = vCandidate;
            }
        }
        m_vPrevStep = v0;
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Force computation time: " << forceTime << "ms\n";
        std::cout << "  Hessian computation time: " << hessianTime << "ms\n";
        std::cout << "  Energy computation time: " << energyTime << "ms\n";
        std::cout << "  Total step time: " << timer.getMilliseconds() << "ms\n\n";
    }
}

void PGD::reset()
{
    BaseSolver::reset();
    m_alpha = s_maxAlpha;
    m_vPrevStep.resize(0);
}

void PGD::getSystemState(VectorX &pos, VectorX &vel) const
{
    const int nObjects = static_cast<int>(m_simObjects.size());
    int accumulatedNDOF = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        auto &object = m_simObjects[objId];
        const int nparticles = object.nParticles();
        for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
            const int startOffset = particleIdx * 3 + accumulatedNDOF;
            pos.segment<3>(startOffset) = object.positions()[particleIdx];
            vel.segment<3>(startOffset) = object.velocities()[particleIdx];
        }
        accumulatedNDOF += object.nDOF();
    }
}
void PGD::setObjectsState(const VectorX &pos, const VectorX &vel, bool updatePos, bool updateVel)
{
    const int nObjects = static_cast<int>(m_simObjects.size());
    int accumulatedNDOF = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        auto &object = m_simObjects[objId];
        const int nparticles = object.nParticles();
        for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
            const int startOffset = particleIdx * 3 + accumulatedNDOF;
            if (updatePos) {
                object.positions()[particleIdx] = pos.segment<3>(startOffset);
            }
            if (updateVel) {
                object.velocities()[particleIdx] = vel.segment<3>(startOffset);
            }
        }
        accumulatedNDOF += object.nDOF();
    }
}
void PGD::getSystemForce(VectorX &f) const
{
    f.setZero();
    const int nObjects = static_cast<int>(m_simObjects.size());
    int accumulatedNDOF = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        const auto &object = m_simObjects[objId];
        // Accumulate internal energy forces
        const auto &energies = object.energies();
        for (auto &energy : energies) {
            const auto nstencils = energy->nStencils();
#pragma omp parallel for
            for (int i = 0; i < nstencils; ++i) {
                energy->accumulateForces(i, object, accumulatedNDOF, f, true);
            }
        }
        accumulatedNDOF += object.nDOF();
    }
}

void PGD::getSystemDiagStiffness(VectorX &K) const
{
    K.setZero();
    const int nObjects = static_cast<int>(m_simObjects.size());
    int accumulatedNDOF = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        const auto &object = m_simObjects[objId];
        // Accumulate internal energy hessians
        const auto &energies = object.energies();
        for (auto &energy : energies) {
            const auto nstencils = energy->nStencils();
#pragma omp parallel for
            for (int i = 0; i < nstencils; ++i) {
                energy->accumulateNegativeHessianDiagonalAtomic(i, object, accumulatedNDOF, K);
            }
        }
        accumulatedNDOF += object.nDOF();
    }
}

void PGD::getSystemDiagMass(VectorX &M) const
{
    const int nObjects = static_cast<int>(m_simObjects.size());
    int accumulatedNDOF = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        auto &object = m_simObjects[objId];
        const auto &masses = object.masses();
        const int nParticles = object.nParticles();
#pragma omp parallel for
        for (int i = 0; i < nParticles; ++i) {
            M.segment<3>(i * 3 + accumulatedNDOF).array() = masses[i];
        }
        accumulatedNDOF += object.nDOF();
    }
}

Real PGD::getSystemEnergy(const VectorX &v, const VectorX &vTilde, const VectorX &M)
{
    Real e = 0;
    const int nObjects = static_cast<int>(m_simObjects.size());
    for (int objId = 0; objId < nObjects; ++objId) {
        auto &object = m_simObjects[objId];
        const auto &energies = object.energies();
        for (auto &energy : energies) {
            const auto nstencils = energy->nStencils();
#pragma omp parallel
            {
                Real threadLocalEnergy = 0;
#pragma omp for
                for (int i = 0; i < nstencils; ++i) {
                    threadLocalEnergy += energy->energy(i, object);
                }
#pragma omp atomic
                e += threadLocalEnergy;
            }
        }
    }
    const VectorX vTerm = v - vTilde;
    e += 0.5 * (vTerm).dot(M.cwiseProduct(vTerm));
    return e;
}
}  // namespace spg::solver