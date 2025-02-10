#pragma once

#include <spg/sim/solver/baseSolver.h>

#include <vector>

namespace spg
{
class SimObject;

namespace solver
{
// Parallel CPU version of the VBD solver by Chen 2024
// (ref: "Vertex Block Descent")
class VBD : public BaseSolver
{
    using ElementEntry = Int2;  // An entry contains the stencil Idx and vertex Idx inside that stencil
    using ElementEntries = std::vector<ElementEntry>;
    using ElementEntriesPerEnergy = std::vector<ElementEntries>;
    using ElementsPerVertex = std::vector<ElementEntriesPerEnergy>;

public:
    enum class InitialGuessType { Inertial, InertialWithAcceleration, Adaptive };
    virtual void step() override;
    virtual void reset() override;
    void setInitialGuessType(const InitialGuessType type) { m_initialGuessType = type; }
    void requirePrecomputationUpdate();

protected:
    void computeParallelVertexGroups();
    void computeStencilInfoPerVertex();

    InitialGuessType m_initialGuessType{InitialGuessType::InertialWithAcceleration};
    int m_iterations{5};  // TODO: Expose this in the solver API
    std::vector<std::vector<std::vector<int>>> m_simObjectsVertexGroups;
    std::vector<ElementsPerVertex> m_simObjectsElementsPerVertex;
    bool m_precomputationUpdateRequired{true};

    std::vector<std::vector<Vector3>> m_simObjPrevStepVelocities;
    // Member structures to avoid unrequired memory reallocations
    std::vector<std::vector<Vector3>> m_simObjectsOldPos;
    std::vector<std::vector<Vector3>> m_simObjectsInertialPositions;
};
}  // namespace solver
}  // namespace spg