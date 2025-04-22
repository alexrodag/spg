#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{

class RBSpringSquaredEnergy : public StencilBlockEnergy<2, 1, RigidBodyGroup>
{
public:
    RBSpringSquaredEnergy() { m_name = "RBSpringSquaredEnergy"; }
    void addStencil(std::array<int, s_stencilSize> stencil,
                    const std::array<Vector3, 2> &localRigidBodyPoints,
                    Real restLength,
                    Real stiffness);
    const std::vector<std::array<Vector3, 2>> &localRBPoints() const { return m_localRBPoints; }
    const std::vector<Real> &squaredRestLengths() const { return m_squaredRestLength; }

protected:
    virtual void dConstraints(int i, const RigidBodyGroup &rbGroup, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const RigidBodyGroup &rbGroup, ConstraintsAD2 &dC) const;
    std::vector<Real> m_squaredRestLength;
    std::vector<std::array<Vector3, 2>> m_localRBPoints;
};
}  // namespace spg