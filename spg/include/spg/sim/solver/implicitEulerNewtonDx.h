#pragma once

#include <spg/sim/solver/implicitEulerBase.h>

#include <vector>
namespace spg
{
namespace solver
{
// Equivalent to ImplicitEulerNewtonDv but using positions as main variable instead of velocities. Implemented just to
// validate that both approaches are numerically equivalent
class ImplicitEulerNewtonDx : public ImplicitEulerBase
{
public:
    virtual void step();
};
}  // namespace solver
}  // namespace spg