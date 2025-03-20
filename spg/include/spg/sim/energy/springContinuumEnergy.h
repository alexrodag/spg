#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
// Uniaxial elastic stretch strain energy as in Spillman 2007
// (ref: "CORDE: Cosserat Rod Elements for the Dynamic Simulation of One-Dimensional Elastic pGroupects")
class SpringContinuumEnergy : public StencilBlockEnergy<2>
{
public:
    SpringContinuumEnergy() { m_name = "SpringContinuumEnergy"; }
    void addStencil(const std::array<int, s_stencilSize> &stencil, const Real restLength, const Real stiffness);
    const std::vector<Real> &restLengths() const { return m_restLength; }

protected:
    virtual void dConstraints(int i, const ParticleGroup &pGroup, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const ParticleGroup &pGroup, ConstraintsAD2 &dC) const;
    std::vector<Real> m_restLength;
};
}  // namespace spg