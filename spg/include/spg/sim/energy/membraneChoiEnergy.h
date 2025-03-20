#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
// Membrane energy corresponding to Baraff 1998, but with the shear component as in Choi 2003
// (refs: "Large Steps in Cloth Simulation", "Extending the Immediate Buckling Model to Triangular Meshes for Simulating
// Complex Clothes")
// The core difference is that it makes shear also linear in x, thus avoiding a mismatch in linearity behabior between
// stretch and shear
class MembraneChoiEnergy : public StencilBlockEnergy<3, 4>
{
public:
    MembraneChoiEnergy() { m_name = "MembraneChoiEnergy"; }
    void addStencil(const std::array<int, s_stencilSize> &stencil, Real kWeft, Real kWarp, Real kBias);
    const std::vector<Matrix2> &inverseReferenceMats() const { return m_inverseReferenceMat; }
    virtual void preparePrecomputations(const ParticleGroup &obj);

protected:
    virtual void dConstraints(int i, const ParticleGroup &obj, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const ParticleGroup &obj, ConstraintsAD2 &dC) const;
    std::vector<Matrix2> m_inverseReferenceMat;
};

}  // namespace spg