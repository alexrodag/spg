#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
class SpringEnergy : public StencilBlockEnergy<2>
{
public:
    SpringEnergy() { m_name = "SpringEnergy"; }
    // Example of manually implemented solver-specific method
    virtual void projectPosition(int i, SimObject &obj, Real dt) const;
    void addStencil(const std::array<int, s_stencilSize> &stencil, const Real restLength, const Real stiffness);
    const std::vector<Real> &restLengths() const { return m_restLength; }

protected:
    // Example of a manually implemented energy
    virtual Real energy(int i, const SimObject &obj) const;
    // Example of a manually implemented energy gradient
    virtual EnergyGradType energyGradient(int i, const SimObject &obj) const;
    // Example of a manually implemented energy hessian
    virtual EnergyHessType energyHessian(int i, const SimObject &obj) const;
    virtual void dEnergy(int i, const SimObject &obj, DScalarFirstD &dE) const;
    virtual void dEnergy(int i, const SimObject &obj, DScalarSecondD &dE) const;
    virtual void dConstraints(int i, const SimObject &obj, DConstraintsFirstD &dC) const;
    virtual void dConstraints(int i, const SimObject &obj, DConstraintsSecondD &dC) const;
    std::vector<Real> m_restLength;
};
}  // namespace spg