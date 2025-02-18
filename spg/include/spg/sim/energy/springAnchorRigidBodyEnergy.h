#pragma once

#include <spg/sim/energy/energy.h>
#include <spg/sim/rigidBodyGroup.h>

namespace spg
{

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
    virtual Real energy(int i, const RigidBodyGroup &obj) const;
    virtual EnergyGrad energyGradient(int i, const RigidBodyGroup &obj) const;
    virtual EnergyHess energyHessian(int i, const RigidBodyGroup &obj) const;
    std::vector<Vector3> m_anchor;
    std::vector<Vector3> m_localRBPoint;
};
}  // namespace spg