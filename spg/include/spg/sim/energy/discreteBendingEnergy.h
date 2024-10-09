#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
// Bending energy corresponding to Rasmus Tamstorf 2013 paper
// (ref: "Derivation of discrete bending forces and their gradients")
class DiscreteBendingEnergy : public StencilBlockEnergy<4>
{
public:
    DiscreteBendingEnergy() { m_name = "DiscreteBendingEnergy"; }
    void addStencil(const std::array<int, s_stencilSize> &stencil, Real restTheta, Real stiffness);
    virtual void preparePrecomputations(const SimObject &obj);
    const std::vector<Real> &restThetas() const { return m_restTheta; }

protected:
    virtual void dConstraints(int i, const SimObject &obj, DConstraintsFirstD &dC) const;
    virtual void dConstraints(int i, const SimObject &obj, DConstraintsSecondD &dC) const;
    std::vector<Real> m_restTheta;
};

}  // namespace spg