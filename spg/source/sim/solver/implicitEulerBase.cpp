#include <spg/sim/solver/implicitEulerBase.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/utils/timer.h>
#include <tbb/tbb.h>

#include <iostream>

#define USE_TBB

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
}

void ImplicitEulerBase::getSystemState(VectorX &pos, VectorX &vel) const
{
    getSystemPositions(pos);
    getSystemVelocities(vel);
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
}
void ImplicitEulerBase::setObjectsState(const VectorX &pos, const VectorX &vel)
{
    setObjectsPositions(pos);
    setObjectsVelocities(vel);
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
#ifdef USE_TBB
            tbb::parallel_for(tbb::blocked_range<int>(0, nstencils),
                              [&energy, &object, &accumulatedNDOF, &f](const tbb::blocked_range<int> &r) {
                                  for (int i = r.begin(); i != r.end(); ++i) {
                                      energy->accumulateForces(i, object, accumulatedNDOF, f, true);
                                  }
                              });
#else
#pragma omp parallel for
            for (int i = 0; i < nstencils; ++i) {
                energy->accumulateForces(i, object, accumulatedNDOF, f, true);
            }
#endif
            timer.stop();
            accumulatedTime += timer.getMilliseconds();
        }
        accumulatedNDOF += object.nDOF();
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
#ifdef USE_TBB
            tbb::enumerable_thread_specific<std::vector<Triplet>> localTriplets;
            tbb::parallel_for(tbb::blocked_range<int>(0, nstencils),
                              [&energy, &object, &accumulatedNDOF, &localTriplets](const tbb::blocked_range<int> &r) {
                                  for (int i = r.begin(); i != r.end(); ++i) {
                                      energy->appendHessianTriplets(i, object, accumulatedNDOF, localTriplets.local());
                                  }
                              });
#else
#pragma omp parallel for
            for (int i = 0; i < nstencils; ++i) {
                perStencilTriplets[i] = energy->negativeHessianTriplets(i, object, accumulatedNDOF);
            }
#endif
            timer.stop();
            accumulatedTimeHComp += timer.getMilliseconds();
            timer.start();
#ifdef USE_TBB
            localTriplets.combine_each(
                [this](const auto &x) { m_tripletHolder.insert(m_tripletHolder.end(), std::begin(x), std::end(x)); });
#else
            for (int i = 0; i < nstencils; ++i) {
                m_tripletHolder.insert(
                    m_tripletHolder.end(), perStencilTriplets[i].begin(), perStencilTriplets[i].end());
            }
#endif
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