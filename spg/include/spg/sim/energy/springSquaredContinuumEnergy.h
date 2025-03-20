#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
// Squared uniaxial continuum spring potential
class SpringSquaredContinuumEnergy : public StencilBlockEnergy<2>
{
public:
    SpringSquaredContinuumEnergy() { m_name = "SpringSquaredContinuumEnergy"; }
    void addStencil(const std::array<int, s_stencilSize> &stencil, const Real restLength, const Real stiffness);
    const std::vector<Real> &restLengths() const { return m_restLength; }

protected:
    virtual void dConstraints(int i, const ParticleGroup &pGroup, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const ParticleGroup &pGroup, ConstraintsAD2 &dC) const;
    std::vector<Real> m_restLength;
};
}  // namespace spg