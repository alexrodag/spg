#pragma once

#include <spg/sim/solver/baseSolver.h>

#include <vector>

namespace spg
{
namespace solver
{

// Experimental, custom Jacobi-preconditioned gradient descent solver, following some ideas from
// refs: "Descent Methods for Elastic Body Simulation on the GPU" and
// "Primal/Dual Descent Methods for Dynamics"
class PGD : public BaseSolver
{
public:
    enum class InitialGuessType { Inertial, InertialWithAcceleration, Adaptive };
    virtual void step() override;
    virtual void reset() override;
    void setInitialGuessType(const InitialGuessType type) { m_initialGuessType = type; }
    void setIterations(const int iterations) { m_iterations = iterations; }

protected:
    void getSystemState(VectorX &pos, VectorX &vel) const;
    void setObjectsState(const VectorX &pos, const VectorX &vel, bool updatePos = true, bool updateVel = true);
    void getSystemForce(VectorX &f) const;
    void getSystemDiagStiffness(VectorX &K) const;
    void getSystemDiagMass(VectorX &M) const;
    Real getSystemEnergy(const VectorX &v, const VectorX &vi, const VectorX &M);

protected:
    static constexpr Real s_maxAlpha{0.5};
    int m_preconditionerUpdateCadency{30};
    int m_iterations{60};  // TODO: Expose this in the solver API
    InitialGuessType m_initialGuessType{InitialGuessType::InertialWithAcceleration};
    Real m_alpha{s_maxAlpha};
    Real m_alphaSplittingFactor{0.75};
    VectorX m_preconditioner;
    VectorX m_vPrevStep;
};
}  // namespace solver
}  // namespace spg