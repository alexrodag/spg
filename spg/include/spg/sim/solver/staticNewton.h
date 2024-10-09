#pragma once

#include <spg/sim/solver/implicitEulerBase.h>

#include <vector>
namespace spg
{
namespace solver
{
// Implicit static solver for the F = 0 problem
class StaticNewton : public ImplicitEulerBase
{
public:
    virtual void step();
};
}  // namespace solver
}  // namespace spg