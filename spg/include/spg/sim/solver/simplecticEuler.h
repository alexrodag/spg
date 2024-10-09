#pragma once

#include <spg/sim/solver/baseSolver.h>

#include <vector>
namespace spg
{
namespace solver
{
// Explicit solver in its simplectic form. Requires really low time step (or equivalently high number of substeps) to
// avoid divergence
class SimplecticEuler : public BaseSolver
{
public:
    virtual void step() override;
};
}  // namespace solver
}  // namespace spg