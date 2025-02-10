#pragma once

#include <spg/sim/solver/baseSolver.h>

#include <vector>
namespace spg
{
template <class TSimObject>
class Energy;
class SimObject;

namespace solver
{
// Parallel CPU XPBD solver (ref: "XPBD: Position-Based Simulation of Compliant Constrained Dynamics") using the single
// iteration per substep strategy (ref: "Small Steps in Physics Simulation")
class XPBD : public BaseSolver
{
public:
    XPBD(bool useParallelGaussSeidel = true);
    virtual void step() override;

protected:
    void computeParallelStencilGroups();

    std::vector<std::vector<Vector3>> m_simObjectsOldPos;
    // TODO: How could this be a more general container of all possible Energy specializations
    std::unordered_map<Energy<SimObject> *, std::vector<std::vector<int>>> m_stencilGroupsPerEnergy;
    bool m_useParallelGaussSeidel;
    bool m_parallelGroupsDirty{true};
};
}  // namespace solver
}  // namespace spg