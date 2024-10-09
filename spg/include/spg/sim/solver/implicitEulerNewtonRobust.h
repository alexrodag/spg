#pragma once

#include <spg/sim/solver/implicitEulerBase.h>

#include <vector>
namespace spg
{
namespace solver
{
// Robust implicit solver based on the ideas presented in
// ref: "Seeking Efficiency for the Accurate Draping of Digital Garments in Production", 2024
class ImplicitEulerNewtonRobust : public ImplicitEulerBase
{
public:
    virtual void step();
    virtual void reset();

protected:
    const Real m_newtonThreshold{1e-5};

    Real m_currentDt{-1.};
    int m_consecutiveSuccessfulSteps{0};
    const int m_requiredSuccessfulStepsToIncreaseDt{5};
    const Real m_timeSplittingFactor{0.75};
    const int m_maxNewtonIterations{5};
    const int m_maxLineSearchIterations{5};
};
}  // namespace solver
}  // namespace spg