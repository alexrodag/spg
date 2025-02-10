#include <spg/sim/solver/implicitEulerBase.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/utils/timer.h>

#include <iostream>

namespace spg::solver
{
void ImplicitEulerBase::getSystemPositions(VectorX &pos) const
{
    const int nObjects = static_cast<int>(m_simObjects.size());
    int accumulatedNDOF = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        const auto &object = m_simObjects[objId];
        const int nparticles = object.nParticles();
        for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
            const int startOffset = particleIdx * 3 + accumulatedNDOF;
            pos.segment<3>(startOffset) = object.positions()[particleIdx];
        }
        accumulatedNDOF += object.nDOF();
    }
    const int nRBs = static_cast<int>(m_rigidBodyGroups.size());
    for (int objId = 0; objId < nRBs; ++objId) {
        const auto &object = m_rigidBodyGroups[objId];
        const int nbodies = object.nBodies();
        for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
            const int startOffset = bodyIdx * 6 + accumulatedNDOF;
            pos.segment<3>(startOffset) = object.positions()[bodyIdx];
            pos.segment<3>(startOffset + 3) = object.orientations()[bodyIdx];
        }
        accumulatedNDOF += object.nDOF();
    }
}

void ImplicitEulerBase::getSystemVelocities(VectorX &vel) const
{
    const int nObjects = static_cast<int>(m_simObjects.size());
    int accumulatedNDOF = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        const auto &object = m_simObjects[objId];
        const int nparticles = object.nParticles();
        for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
            const int startOffset = particleIdx * 3 + accumulatedNDOF;
            vel.segment<3>(startOffset) = object.velocities()[particleIdx];
        }
        accumulatedNDOF += object.nDOF();
    }
    const int nRBs = static_cast<int>(m_rigidBodyGroups.size());
    for (int objId = 0; objId < nRBs; ++objId) {
        const auto &object = m_rigidBodyGroups[objId];
        const int nbodies = object.nBodies();
        for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
            const int startOffset = bodyIdx * 6 + accumulatedNDOF;
            vel.segment<3>(startOffset) = object.velocities()[bodyIdx];
            vel.segment<3>(startOffset + 3) = object.omegas()[bodyIdx];
        }
        accumulatedNDOF += object.nDOF();
    }
}

void ImplicitEulerBase::setObjectsPositions(const VectorX &pos)
{
    const int nObjects = static_cast<int>(m_simObjects.size());
    int accumulatedNDOF = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        auto &object = m_simObjects[objId];
        const int nparticles = object.nParticles();
        for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
            const int startOffset = particleIdx * 3 + accumulatedNDOF;
            object.positions()[particleIdx] = pos.segment<3>(startOffset);
        }
        accumulatedNDOF += object.nDOF();
    }
    const int nRBs = static_cast<int>(m_rigidBodyGroups.size());
    for (int objId = 0; objId < nRBs; ++objId) {
        auto &object = m_rigidBodyGroups[objId];
        const int nbodies = object.nBodies();
        for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
            const int startOffset = bodyIdx * 6 + accumulatedNDOF;
            object.positions()[bodyIdx] = pos.segment<3>(startOffset);
            object.orientations()[bodyIdx] = pos.segment<3>(startOffset + 3);
        }
        object.updateRotationMatrices();
        object.updateInertias();
        accumulatedNDOF += object.nDOF();
    }
}

void ImplicitEulerBase::setObjectsVelocities(const VectorX &vel)
{
    const int nObjects = static_cast<int>(m_simObjects.size());
    int accumulatedNDOF = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        auto &object = m_simObjects[objId];
        const int nparticles = object.nParticles();
        for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
            const int startOffset = particleIdx * 3 + accumulatedNDOF;
            object.velocities()[particleIdx] = vel.segment<3>(startOffset);
        }
        accumulatedNDOF += object.nDOF();
    }
    const int nRBs = static_cast<int>(m_rigidBodyGroups.size());
    for (int objId = 0; objId < nRBs; ++objId) {
        auto &object = m_rigidBodyGroups[objId];
        const int nbodies = object.nBodies();
        for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
            const int startOffset = bodyIdx * 6 + accumulatedNDOF;
            object.velocities()[bodyIdx] = vel.segment<3>(startOffset);
            object.omegas()[bodyIdx] = vel.segment<3>(startOffset + 3);
        }
        accumulatedNDOF += object.nDOF();
    }
}

void ImplicitEulerBase::updateObjectsStateDv(const VectorX &dv, const Real dt)
{
    // TODO: Move this to the simObjects
    const int nObjects = static_cast<int>(m_simObjects.size());
    int accumulatedNDOF = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        auto &object = m_simObjects[objId];
        const int nparticles = object.nParticles();
        for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
            const int startOffset = particleIdx * 3 + accumulatedNDOF;
            object.velocities()[particleIdx] += dv.segment<3>(startOffset);
            object.positions()[particleIdx] += object.velocities()[particleIdx] * dt;
        }
        accumulatedNDOF += object.nDOF();
    }
    const int nRBs = static_cast<int>(m_rigidBodyGroups.size());
    for (int objId = 0; objId < nRBs; ++objId) {
        auto &object = m_rigidBodyGroups[objId];
        const int nbodies = object.nBodies();
        for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
            const int startOffset = bodyIdx * 6 + accumulatedNDOF;
            // linear part
            object.velocities()[bodyIdx] += dv.segment<3>(startOffset);
            object.positions()[bodyIdx] += object.velocities()[bodyIdx] * dt;
            // angular part
            Vector3 &omega = object.omegas()[bodyIdx];
            Matrix3 &rotMatrix = object.rotationMatrices()[bodyIdx];
            omega += dv.segment<3>(startOffset + 3);
            if (const auto omegaNorm = omega.norm(); omegaNorm != 0) {
                // Incremental compositions of the rotation
                rotMatrix =
                    Eigen::AngleAxis<spg::Real>(omegaNorm * dt, omega / omegaNorm).toRotationMatrix() * rotMatrix;
            }
        }
        // Store the current orientation in the thetas and update inertias
        object.updateThetas();
        object.updateInertias();
        accumulatedNDOF += object.nDOF();
    }
}

void ImplicitEulerBase::updateObjectsStateDx(const VectorX &dx, const Real dt)
{
}

void ImplicitEulerBase::getSystemForce(VectorX &f) const
{
    f.setZero();
    const int nObjects = static_cast<int>(m_simObjects.size());
    int accumulatedNDOF = 0;
    Timer timer;
    Real accumulatedTime = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        const auto &object = m_simObjects[objId];
        // Accumulate external forces (gravity)
        const auto &masses = object.masses();
        const int nMasses = static_cast<int>(masses.size());
        for (int i = 0; i < nMasses; ++i) {
            const Real mass = masses[i];
            f.segment<3>(i * 3 + accumulatedNDOF) += m_gravity * mass;
        }
        // Accumulate internal energy forces
        // TODO: Forces are being computed two times in some cases, try to use the gradients from the hessian
        // computation to avoid it
        const auto &energies = object.energies();
        for (const auto &energy : energies) {
            const auto nstencils = energy->nStencils();
            timer.start();
#pragma omp parallel for
            for (int i = 0; i < nstencils; ++i) {
                energy->accumulateForces(i, object, accumulatedNDOF, f, true);
            }
            timer.stop();
            accumulatedTime += timer.getMilliseconds();
        }
        accumulatedNDOF += object.nDOF();
    }

    const int nRigidBodyGroups = static_cast<int>(m_rigidBodyGroups.size());
    for (int objId = 0; objId < nRigidBodyGroups; ++objId) {
        const auto &rigidBodyGroup = m_rigidBodyGroups[objId];
        // Accumulate external torques and fictional torques
        const auto &inertias = rigidBodyGroup.inertias();
        const int nbodies = static_cast<int>(rigidBodyGroup.nBodies());

        auto l_skew = [](const spg::Vector3 &v) {
            spg::Matrix3 vSkew;
            vSkew << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
            return vSkew;
        };
        // Add gravity forces and gyro torques
        for (int i = 0; i < nbodies; ++i) {
            f.segment<3>(i * 6 + accumulatedNDOF) += m_gravity * rigidBodyGroup.masses()[i];
            f.segment<3>(i * 6 + 3 + accumulatedNDOF) +=
                -l_skew(rigidBodyGroup.omegas()[i]) * rigidBodyGroup.inertias()[i] * rigidBodyGroup.omegas()[i];
        }
        // Accumulate internal energy forces
        // TODO: Forces are being computed two times in some cases, try to use the gradients from the hessian
        // computation to avoid it
        const auto &energies = rigidBodyGroup.energies();
        for (const auto &energy : energies) {
            const auto nstencils = energy->nStencils();
            timer.start();
#pragma omp parallel for
            for (int i = 0; i < nstencils; ++i) {
                energy->accumulateForces(i, rigidBodyGroup, accumulatedNDOF, f, true);
            }
            timer.stop();
            accumulatedTime += timer.getMilliseconds();
        }
        accumulatedNDOF += rigidBodyGroup.nDOF();
    }

    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  f computation time: " << accumulatedTime << "ms\n";
    }
}

void ImplicitEulerBase::getSystemMassMatrix(SparseMatrix &M) const
{
    Timer timer;
    timer.start();
    const int nObjects = static_cast<int>(m_simObjects.size());
    m_tripletHolder.clear();
    int accumulatedNDOF = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        const auto &object = m_simObjects[objId];
        // Accumulate mass triplets
        const auto &masses = object.masses();
        for (int i = 0; i < masses.size(); ++i) {
            const Real mass = masses[i];
            for (int j = 0; j < 3; ++j) {
                const int diagIndex = i * 3 + j + accumulatedNDOF;
                m_tripletHolder.emplace_back(diagIndex, diagIndex, mass);
            }
        }

        accumulatedNDOF += object.nDOF();
    }
    const int nRigidBodyGroups = static_cast<int>(m_rigidBodyGroups.size());
    for (int objId = 0; objId < nRigidBodyGroups; ++objId) {
        const auto &object = m_rigidBodyGroups[objId];
        // Accumulate mass triplets
        const auto &masses = object.masses();
        const auto &inertias = object.inertias();
        for (int i = 0; i < masses.size(); ++i) {
            const Real mass = masses[i];
            for (int j = 0; j < 3; ++j) {
                const int diagIndex = i * 6 + j + accumulatedNDOF;
                m_tripletHolder.emplace_back(diagIndex, diagIndex, mass);
            }
            const int upperLeftCornerIndex = i * 6 + 3 + accumulatedNDOF;
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    m_tripletHolder.emplace_back(upperLeftCornerIndex + j, upperLeftCornerIndex + k, inertias[i](j, k));
                }
            }
        }
        accumulatedNDOF += object.nDOF();
    }
    M.setFromTriplets(m_tripletHolder.begin(), m_tripletHolder.end());
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  M computation and assembly time: " << timer.getMilliseconds() << "ms\n";
    }
}

void ImplicitEulerBase::getSystemStiffnessMatrix(SparseMatrix &K) const
{
    const int nObjects = static_cast<int>(m_simObjects.size());
    m_tripletHolder.clear();
    int accumulatedNDOF = 0;
    Timer timer;
    Real accumulatedTimeHComp = 0;
    Real accumulatedTimeHTriplets = 0;
    for (int objId = 0; objId < nObjects; ++objId) {
        const auto &object = m_simObjects[objId];
        // Accumulate hessian triplets
        const auto &energies = object.energies();
        for (const auto &energy : energies) {
            const auto nstencils = energy->nStencils();
            std::vector<std::vector<Triplet>> perStencilTriplets(nstencils);
            timer.start();
#pragma omp parallel for
            for (int i = 0; i < nstencils; ++i) {
                perStencilTriplets[i] = energy->negativeHessianTriplets(i, object, accumulatedNDOF);
            }
            timer.stop();
            accumulatedTimeHComp += timer.getMilliseconds();
            timer.start();
            for (int i = 0; i < nstencils; ++i) {
                m_tripletHolder.insert(
                    m_tripletHolder.end(), perStencilTriplets[i].begin(), perStencilTriplets[i].end());
            }
            timer.stop();
            accumulatedTimeHTriplets += timer.getMilliseconds();
        }
        accumulatedNDOF += object.nDOF();
    }

    const int nRigidBodyGroups = static_cast<int>(m_rigidBodyGroups.size());
    for (int objId = 0; objId < nRigidBodyGroups; ++objId) {
        const auto &object = m_rigidBodyGroups[objId];
        // Accumulate hessian triplets
        const auto &energies = object.energies();
        for (const auto &energy : energies) {
            const auto nstencils = energy->nStencils();
            std::vector<std::vector<Triplet>> perStencilTriplets(nstencils);
            timer.start();
#pragma omp parallel for
            for (int i = 0; i < nstencils; ++i) {
                perStencilTriplets[i] = energy->negativeHessianTriplets(i, object, accumulatedNDOF);
            }
            timer.stop();
            accumulatedTimeHComp += timer.getMilliseconds();
            timer.start();
            for (int i = 0; i < nstencils; ++i) {
                m_tripletHolder.insert(
                    m_tripletHolder.end(), perStencilTriplets[i].begin(), perStencilTriplets[i].end());
            }
            timer.stop();
            accumulatedTimeHTriplets += timer.getMilliseconds();
        }
        accumulatedNDOF += object.nDOF();
    }

    timer.start();
    K.setFromTriplets(m_tripletHolder.begin(), m_tripletHolder.end());
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  H computation time: " << accumulatedTimeHComp << "ms\n";
        std::cout << "  H tripplets accumulation time: " << accumulatedTimeHTriplets << "ms\n";
        std::cout << "  H assembly time: " << timer.getMilliseconds() << "ms\n";
    }
}

void ImplicitEulerBase::solveLinearSystem(const SparseMatrix &A, const VectorX &b, VectorX &x) const
{
    // Note: Under certain circumstances, some direct solver such as SimplicialLLT may be faster
    Timer timer;
    timer.start();
    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> cg;
    cg.setTolerance(m_cgThreshold);
    cg.setMaxIterations(m_cgMaxIterations);
    cg.compute(A);
    x = cg.solve(b);
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Linear solve time: " << timer.getMilliseconds() << "ms. CG iterations: " << cg.iterations()
                  << "\n";
    }
}
}  // namespace spg::solver