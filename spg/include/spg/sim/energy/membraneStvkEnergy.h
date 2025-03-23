#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
// Membrane energy corresponding to Volino 2009
// (ref: "A Simple Approach to Nonlinear Tensile Stiffness for Accurate Cloth Simulation")
class MembraneStvkEnergy : public StencilBlockEnergy<3, 3>
{
public:
    MembraneStvkEnergy() { m_name = "MembraneStvkEnergy"; }
    void addStencil(const std::array<int, s_stencilSize> &stencil, Real young, Real poisson);
    const std::vector<Matrix2> &inverseReferenceMats() const { return m_inverseReferenceMat; }
    virtual void preparePrecomputations(const ParticleGroup &obj);

protected:
    virtual void dConstraints(int i, const ParticleGroup &obj, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const ParticleGroup &obj, ConstraintsAD2 &dC) const;
    std::vector<Matrix2> m_inverseReferenceMat;
};

}  // namespace spg