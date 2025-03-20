#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
// Classic spring energy
// (ref: "https://en.wikipedia.org/wiki/Elastic_energy")
class SpringEnergy : public StencilBlockEnergy<2>
{
public:
    SpringEnergy() { m_name = "SpringEnergy"; }
    // Example of manually implemented solver-specific method
    virtual void projectPosition(int i, ParticleGroup &obj, Real dt) const;
    void addStencil(const std::array<int, s_stencilSize> &stencil, const Real restLength, const Real stiffness);
    const std::vector<Real> &restLengths() const { return m_restLength; }

protected:
    // Example of a manually implemented energy
    virtual Real energy(int i, const ParticleGroup &obj) const;
    // Example of a manually implemented energy gradient
    virtual EnergyGrad energyGradient(int i, const ParticleGroup &obj) const;
    // Example of a manually implemented energy hessian
    virtual EnergyHess energyHessian(int i, const ParticleGroup &obj) const;
    virtual void dEnergy(int i, const ParticleGroup &obj, RealAD1 &dE) const;
    virtual void dEnergy(int i, const ParticleGroup &obj, RealAD2 &dE) const;
    virtual void dConstraints(int i, const ParticleGroup &obj, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const ParticleGroup &obj, ConstraintsAD2 &dC) const;
    std::vector<Real> m_restLength;
};
}  // namespace spg