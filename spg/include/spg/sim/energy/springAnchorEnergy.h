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
    virtual void dEnergy(int i, const SimObject &obj, RealAD1 &dE) const;
    virtual void dEnergy(int i, const SimObject &obj, RealAD2 &dE) const;
    virtual void dConstraints(int i, const SimObject &obj, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const SimObject &obj, ConstraintsAD2 &dC) const;
    std::vector<Vector3> m_anchor;
};

class SpringAnchorRBEnergy : public StencilBlockEnergy<1, 1, RigidBodyGroup>
{
public:
    SpringAnchorRBEnergy() { m_name = "SpringRBAnchorEnergy"; }
    void addStencil(std::array<int, s_stencilSize> stencil,
                    const Vector3 &anchor,
                    const Vector3 &localRigidBodyPoint,
                    Real stiffness);
    const std::vector<Vector3> &anchors() const { return m_anchor; }
    const std::vector<Vector3> &localRBPoints() const { return m_localRBPoint; }

protected:
    virtual void dEnergy(int i, const RigidBodyGroup &obj, RealAD1 &dE) const;
    virtual void dEnergy(int i, const RigidBodyGroup &obj, RealAD2 &dE) const;
    virtual void dConstraints(int i, const RigidBodyGroup &obj, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const RigidBodyGroup &obj, ConstraintsAD2 &dC) const;
    std::vector<Vector3> m_anchor;
    std::vector<Vector3> m_localRBPoint;
};
}  // namespace spg