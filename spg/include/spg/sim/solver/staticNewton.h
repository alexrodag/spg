#pragma once

#include <spg/sim/solver/implicitEulerBase.h>

#include <vector>
namespace spg
{
namespace solver
{
// Implicit static solver for the F = 0 problem
// (Ref: https://crl.ethz.ch/teaching/shape-modeling-18/lectures/10_PhysicsDeformations.pdf)
class StaticNewton : public ImplicitEulerBase
{
public:
    virtual void step();
};
}  // namespace solver
}  // namespace spg