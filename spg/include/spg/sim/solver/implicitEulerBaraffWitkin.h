#pragma once

#include <spg/sim/solver/implicitEulerBase.h>

#include <vector>
namespace spg
{
namespace solver
{
// ref: "Large Steps in Cloth Simulation", 1998
class ImplicitEulerBaraffWitkin : public ImplicitEulerBase
{
public:
    virtual void step();
};
}  // namespace solver
}  // namespace spg