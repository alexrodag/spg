#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{

class SpringAnchorRBEnergy : public StencilBlockEnergy<1, 1, RigidBodyGroup>
{
public:
    SpringAnchorRBEnergy() { m_name = "SpringRBAnchorEnergy"; }
    void addStencil(std::array<int, s_stencilSize> stencil,
                    const Vector3 &localRigidBodyPoint,
                    const Vector3 &anchor,
                    Real stiffness);
    const std::vector<Vector3> &anchors() const { return m_anchor; }
    const std::vector<Vector3> &localRBPoints() const { return m_localRBPoint; }
    void updateAnchor(int i, const Vector3 &anchor) { m_anchor[i] = anchor; }

protected:
    virtual void dEnergy(int i, const RigidBodyGroup &rbGroup, RealAD1 &dE) const;
    virtual void dEnergy(int i, const RigidBodyGroup &rbGroup, RealAD2 &dE) const;
    virtual void dConstraints(int i, const RigidBodyGroup &rbGroup, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const RigidBodyGroup &rbGroup, ConstraintsAD2 &dC) const;
    virtual Real energy(int i, const RigidBodyGroup &rbGroup) const;
    virtual EnergyGrad energyGradient(int i, const RigidBodyGroup &rbGroup) const;
    virtual EnergyHess energyHessian(int i, const RigidBodyGroup &rbGroup) const;
    std::vector<Vector3> m_anchor;
    std::vector<Vector3> m_localRBPoint;
};
}  // namespace spg