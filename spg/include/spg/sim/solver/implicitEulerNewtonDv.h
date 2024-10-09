#pragma once

#include <spg/sim/solver/implicitEulerBase.h>

#include <vector>
namespace spg
{
namespace solver
{
// Simplified version of ImplicitEulerNewtonRobust to do a single newton iteration. Implemented mainly for comparison
// purposes with the BaraffWitking version
class ImplicitEulerNewtonDv : public ImplicitEulerBase
{
public:
    virtual void step();
};
}  // namespace solver
}  // namespace spg