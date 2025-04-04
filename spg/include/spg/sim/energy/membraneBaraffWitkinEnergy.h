#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
// Membrane energy corresponding to Baraff 1998
// (ref: "Large Steps in Cloth Simulation")
class MembraneBaraffWitkinEnergy : public StencilBlockEnergy<3, 3>
{
public:
    MembraneBaraffWitkinEnergy() { m_name = "MembraneBaraffWitkinEnergy"; }
    void addStencil(const std::array<int, s_stencilSize> &stencil, Real kWeft, Real kWarp, Real kBias);
    const std::vector<Matrix2> &inverseReferenceMats() const { return m_inverseReferenceMat; }
    virtual void preparePrecomputations(const ParticleGroup &pGroup);

protected:
    virtual void dConstraints(int i, const ParticleGroup &pGroup, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const ParticleGroup &pGroup, ConstraintsAD2 &dC) const;
    std::vector<Matrix2> m_inverseReferenceMat;
};

}  // namespace spg