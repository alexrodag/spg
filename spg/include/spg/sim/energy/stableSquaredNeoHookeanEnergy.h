#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
// FEM energy corresponding to a variation of the Stable Neo Hookean energy found in the supplementary material of the
// work "Constraint-based Formulation of Stable Neo-Hookean Materials" by Macklin et al. 2021, squaring the Deviatoric
// term. This strategy is also mentioned in the SIGGRAPH course "Dynamic Deformables: Implementation and Production
// Practicalities" by Kim et al. 2022.
class StableSquaredNeoHookeanEnergy : public StencilBlockEnergy<4, 2>
{
public:
    StableSquaredNeoHookeanEnergy() { m_name = "StableSquaredNeoHookeanEnergy"; }
    void addStencil(const std::array<int, s_stencilSize> &stencil, Real young, Real poisson);
    virtual void preparePrecomputations(const ParticleGroup &obj);
    const std::vector<Matrix3> &inverseReferenceMats() const { return m_inverseReferenceMat; }

protected:
    virtual void dConstraints(int i, const ParticleGroup &obj, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const ParticleGroup &obj, ConstraintsAD2 &dC) const;
    std::vector<Matrix3> m_inverseReferenceMat;
};

}  // namespace spg