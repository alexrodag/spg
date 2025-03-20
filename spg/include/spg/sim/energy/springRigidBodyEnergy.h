#pragma once

#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject/rigidBodyGroup.h>

namespace spg
{

class SpringRBEnergy : public StencilBlockEnergy<2, 1, RigidBodyGroup>
{
public:
    SpringRBEnergy() { m_name = "SpringRBEnergy"; }
    void addStencil(std::array<int, s_stencilSize> stencil,
                    const std::array<Vector3, 2> &localRigidBodyPoints,
                    Real restLength,
                    Real stiffness);
    const std::vector<std::array<Vector3, 2>> &localRBPoints() const { return m_localRBPoints; }
    const std::vector<Real> &restLengths() const { return m_restLength; }

protected:
    virtual void dConstraints(int i, const RigidBodyGroup &obj, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const RigidBodyGroup &obj, ConstraintsAD2 &dC) const;
    std::vector<Real> m_restLength;
    std::vector<std::array<Vector3, 2>> m_localRBPoints;
};
}  // namespace spg