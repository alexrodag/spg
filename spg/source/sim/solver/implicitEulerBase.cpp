#include <spg/sim/solver/implicitEulerBase.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/sim/rigidBodyGroup.h>
#include <spg/utils/timer.h>
#include <spg/utils/functionalUtilities.h>

#include <iostream>

namespace spg::solver
{
void ImplicitEulerBase::getSystemPositions(VectorX &pos) const
{
    int accumulatedNDOF = 0;
    apply_each(
        [&accumulatedNDOF, &pos](const auto &objs) {
            for (const auto &obj : objs) {
                obj.getPositions(pos, accumulatedNDOF);
                accumulatedNDOF += obj.nDOF();
            }
        },
        m_objects);
}

void ImplicitEulerBase::getSystemVelocities(VectorX &vel) const
{
    int accumulatedNDOF = 0;
    apply_each(
        [&accumulatedNDOF, &vel](const auto &objs) {
            for (const auto &obj : objs) {
                obj.getVelocities(vel, accumulatedNDOF);
                accumulatedNDOF += obj.nDOF();
            }
        },
        m_objects);
}

void ImplicitEulerBase::setObjectsPositions(const VectorX &pos)
{
    int accumulatedNDOF = 0;
    apply_each(
        [&accumulatedNDOF, &pos](auto &objs) {
            for (auto &obj : objs) {
                obj.setPositions(pos, accumulatedNDOF);
                accumulatedNDOF += obj.nDOF();
            }
        },
        m_objects);
}

void ImplicitEulerBase::setObjectsVelocities(const VectorX &vel)
{
    int accumulatedNDOF = 0;
    apply_each(
        [&accumulatedNDOF, &vel](auto &objs) {
            for (auto &obj : objs) {
                obj.setVelocities(vel, accumulatedNDOF);
                accumulatedNDOF += obj.nDOF();
            }
        },
        m_objects);
}

void ImplicitEulerBase::integrateObjectsVelocities(const Real dt)
{
    apply_each(
        [dt](auto &objs) {
            for (auto &obj : objs) {
                obj.integrateVelocities(dt);
            }
        },
        m_objects);
}

void ImplicitEulerBase::integrateObjectsVelocitiesFromDx(const VectorX &dx, const VectorX &oldPos, const Real invdt)
{
    int accumulatedNDOF = 0;
    apply_each(
        [&accumulatedNDOF, &dx, &oldPos, invdt](auto &objs) {
            for (auto &obj : objs) {
                obj.updatePositionsFromDx(dx, accumulatedNDOF);
                obj.computeIntegratedVelocities(oldPos, accumulatedNDOF, invdt);
                accumulatedNDOF += obj.nDOF();
            }
        },
        m_objects);
}

void ImplicitEulerBase::updateObjectsPositionsFromDx(const VectorX &dx)
{
    int accumulatedNDOF = 0;
    apply_each(
        [&accumulatedNDOF, &dx](auto &objs) {
            for (auto &obj : objs) {
                obj.updatePositionsFromDx(dx, accumulatedNDOF);
                accumulatedNDOF += obj.nDOF();
            }
        },
        m_objects);
}

void ImplicitEulerBase::getSystemForce(VectorX &f) const
{
    int accumulatedNDOF = 0;
    Timer timer;
    Real accumulatedTime = 0;
    f.setZero();
    apply_each(
        [this, &accumulatedNDOF, &f, &timer, &accumulatedTime](auto &objs) {
            auto l_skew = [](const spg::Vector3 &v) {
                spg::Matrix3 vSkew;
                vSkew << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
                return vSkew;
            };

            for (const auto &obj : objs) {
                const auto &masses = obj.masses();
                const int nPrimitives = static_cast<int>(obj.size());
                for (int i = 0; i < nPrimitives; ++i) {
                    const Real mass = masses[i];
                    f.segment<3>(i * obj.s_nDOFs + accumulatedNDOF) += m_gravity * mass;
                    if constexpr (std::is_same_v<std::decay_t<decltype(obj)>, RigidBodyGroup>) {
                        // TODO: Transform this into an implicit energy, as it can lead to instability
                        f.segment<3>(i * obj.s_nDOFs + 3 + accumulatedNDOF) +=
                            -l_skew(obj.omegas()[i]) * obj.inertias()[i] * obj.omegas()[i];
                    }
                }
                // Accumulate internal energy forces
                // TODO: Forces are being computed two times in some cases, try to use the gradients from the hessian
                // computation to avoid it
                const auto &energies = obj.energies();
                for (const auto &energy : energies) {
                    const auto nstencils = energy->nStencils();
                    timer.start();
                    // #pragma omp parallel for
                    for (int i = 0; i < nstencils; ++i) {
                        energy->accumulateForces(i, obj, accumulatedNDOF, f, true);
                    }
                    timer.stop();
                    accumulatedTime += timer.getMilliseconds();
                }
                accumulatedNDOF += obj.nDOF();
            }
        },
        m_objects);
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  f computation time: " << accumulatedTime << "ms\n";
    }
}

void ImplicitEulerBase::getSystemMassMatrix(SparseMatrix &M) const
{
    Timer timer;
    timer.start();
    int accumulatedNDOF = 0;
    m_tripletHolder.clear();
    apply_each(
        [this, &accumulatedNDOF, &M](auto &objs) {
            for (auto &obj : objs) {
                const int nPrimitives = obj.size();
                for (int i = 0; i < nPrimitives; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        const int diagIndex = i * obj.s_nDOFs + j + accumulatedNDOF;
                        m_tripletHolder.emplace_back(diagIndex, diagIndex, obj.masses()[i]);
                    }
                    if constexpr (std::is_same_v<std::decay_t<decltype(obj)>, RigidBodyGroup>) {
                        const int upperLeftCornerIndex = i * obj.s_nDOFs + 3 + accumulatedNDOF;
                        for (int j = 0; j < 3; ++j) {
                            for (int k = 0; k < 3; ++k) {
                                m_tripletHolder.emplace_back(
                                    upperLeftCornerIndex + j, upperLeftCornerIndex + k, obj.inertias()[i](j, k));
                            }
                        }
                    }
                }
                accumulatedNDOF += obj.nDOF();
            }
        },
        m_objects);
    M.setFromTriplets(m_tripletHolder.begin(), m_tripletHolder.end());
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  M computation and assembly time: " << timer.getMilliseconds() << "ms\n";
    }
}

void ImplicitEulerBase::getSystemStiffnessMatrix(SparseMatrix &K) const
{
    m_tripletHolder.clear();
    int accumulatedNDOF = 0;
    Timer timer;
    Real accumulatedTimeHComp = 0;
    Real accumulatedTimeHTriplets = 0;
    apply_each(
        [this, &accumulatedNDOF, &timer, &accumulatedTimeHComp, &accumulatedTimeHTriplets](auto &objs) {
            for (const auto &obj : objs) {
                // Accumulate hessian triplets
                const auto &energies = obj.energies();
                for (const auto &energy : energies) {
                    const auto nstencils = energy->nStencils();
                    std::vector<std::vector<Triplet>> perStencilTriplets(nstencils);
                    timer.start();
#pragma omp parallel for
                    for (int i = 0; i < nstencils; ++i) {
                        perStencilTriplets[i] = energy->negativeHessianTriplets(i, obj, accumulatedNDOF);
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
                accumulatedNDOF += obj.nDOF();
            }
        },
        m_objects);
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