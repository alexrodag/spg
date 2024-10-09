#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
class SpringAnchorEnergy : public StencilBlockEnergy<1>
{
public:
    SpringAnchorEnergy() { m_name = "SpringAnchorEnergy"; }
    void addStencil(std::array<int, s_stencilSize> stencil, const Vector3 &anchor, Real stiffness);
    const std::vector<Vector3> &anchors() const { return m_anchor; }

protected:
    virtual void dEnergy(int i, const SimObject &obj, DScalarFirstD &dE) const;
    virtual void dEnergy(int i, const SimObject &obj, DScalarSecondD &dE) const;
    virtual void dConstraints(int i, const SimObject &obj, DConstraintsFirstD &dC) const;
    virtual void dConstraints(int i, const SimObject &obj, DConstraintsSecondD &dC) const;
    std::vector<Vector3> m_anchor;
};
}  // namespace spg