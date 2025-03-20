#pragma once

#include <spg/sim/solver/baseSolver.h>

#include <vector>

namespace spg
{
class ParticleGroup;

namespace solver
{
// Parallel CPU version of the VBD solver by Chen 2024
// (ref: "Vertex Block Descent")
// Note: When dealing with rigid bodies, the logic of a "vertex" extends to a rigid body
class VBD : public BaseSolver
{
    using ElementEntry = Int2;  // An entry contains the stencil Idx and vertex Idx inside that stencil
    using ElementEntries = std::vector<ElementEntry>;
    using ElementEntriesPerEnergy = std::vector<ElementEntries>;
    using ElementsPerVertex = std::vector<ElementEntriesPerEnergy>;

public:
    enum class InitialGuessType { InertialWithAcceleration };
    virtual void step() override;
    virtual void reset() override;
    void setInitialGuessType(const InitialGuessType type) { m_initialGuessType = type; }
    void requirePrecomputationUpdate() { m_precomputationUpdateRequired = true; }

protected:
    void computeParallelVertexGroups();
    void computeStencilInfoPerVertex();

    InitialGuessType m_initialGuessType{InitialGuessType::InertialWithAcceleration};
    int m_iterations{1};  // TODO: Expose this in the solver API
    std::vector<std::vector<std::vector<int>>> m_simObjectsVertexGroups;
    std::vector<ElementsPerVertex> m_simObjectsElementsPerVertex;
    bool m_precomputationUpdateRequired{true};

    // Member structures to avoid unrequired memory reallocations
    VectorX m_xOld;
    VectorX m_xInertial;
};
}  // namespace solver
}  // namespace spg