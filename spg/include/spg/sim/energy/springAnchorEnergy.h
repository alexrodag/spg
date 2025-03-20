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
    void updateAnchor(int i, const Vector3 &anchor) { m_anchor[i] = anchor; }

protected:
    virtual void dEnergy(int i, const ParticleGroup &pGroup, RealAD1 &dE) const;
    virtual void dEnergy(int i, const ParticleGroup &pGroup, RealAD2 &dE) const;
    virtual void dConstraints(int i, const ParticleGroup &pGroup, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const ParticleGroup &pGroup, ConstraintsAD2 &dC) const;
    std::vector<Vector3> m_anchor;
};
}  // namespace spg