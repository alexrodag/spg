#pragma once

#include <spg/sim/solver/baseSolver.h>

#include <vector>
namespace spg
{
namespace solver
{
class ImplicitEulerBase : public BaseSolver
{
protected:
    void getSystemPositions(VectorX &pos) const;
    void getSystemVelocities(VectorX &vel) const;
    void setObjectsPositions(const VectorX &pos);
    void setObjectsVelocities(const VectorX &vel);
    void updateObjectsStateDv(const VectorX &dv, const Real dt);
    void updateObjectsStateDx(const VectorX &dx, const Real dt);
    // TODO Evaluate if these methods can use copy elision to improve readability
    void getSystemForce(VectorX &f) const;
    void getSystemStiffnessMatrix(SparseMatrix &K) const;
    void getSystemMassMatrix(SparseMatrix &M) const;
    void solveLinearSystem(const SparseMatrix &A, const VectorX &b, VectorX &x) const;

protected:
    const Real m_cgThreshold{1e-5};
    const int m_cgMaxIterations{300};
    mutable std::vector<Triplet> m_tripletHolder;  // Dynamic holder to avoid reallocation across steps
};
}  // namespace solver
}  // namespace spg