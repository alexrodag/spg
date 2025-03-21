#pragma once

#include <spg/sim/solver/implicitEulerBase.h>

#include <vector>
namespace spg
{
namespace solver
{
// Quasi-static solver with mass-based regularization based on the ideas presented in Zhang 2022
// (ref: "Progressive Simulation for Cloth Quasistatics")
class QuasiStaticNewtonRobust : public ImplicitEulerBase
{
public:
    virtual void step();
    virtual void reset();

protected:
    void computeSystemIntegratedVelocities(const VectorX &x0, VectorX &v, Real invdt);

    const Real m_newtonThreshold{1e-2};

    // Quality status control
    Real m_currentDt{-1.};
    int m_consecutiveSuccessfulSteps{0};
    const int m_requiredSuccessfulStepsToIncreaseDt{5};
    const Real m_timeSplittingFactor{0.75};
    const int m_maxNewtonIterations{5};
    const int m_maxLineSearchIterations{5};
};
}  // namespace solver
}  // namespace spg