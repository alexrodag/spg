#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{

class RBOrientationAlignmentEnergy : public StencilBlockEnergy<2, 1, RigidBodyGroup>
{
public:
    RBOrientationAlignmentEnergy() { m_name = "RBOrientationAlignmentEnergy"; }
    void addStencil(std::array<int, s_stencilSize> stencil, Real stiffness);

protected:
    virtual void dConstraints(int i, const RigidBodyGroup &rbGroup, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const RigidBodyGroup &rbGroup, ConstraintsAD2 &dC) const;
};
}  // namespace spg