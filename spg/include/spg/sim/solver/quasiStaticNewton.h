#pragma once

#include <spg/sim/solver/implicitEulerBase.h>

#include <vector>
namespace spg
{
namespace solver
{
// Simplified version of QuasiStaticNewtonRobust to do a single newton iteration.
class QuasiStaticNewton : public ImplicitEulerBase
{
public:
    virtual void step();
};
}  // namespace solver
}  // namespace spg