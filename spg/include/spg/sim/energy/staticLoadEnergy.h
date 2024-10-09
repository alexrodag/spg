#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
class StaticLoadEnergy : public StencilBlockEnergy<1>
{
public:
    StaticLoadEnergy() { m_name = "StaticLoadEnergy"; }
    void addStencil(std::array<int, s_stencilSize> stencil, const Vector3 &loadDirection, const Real stiffness)
    {
        m_stencils.push_back(stencil);
        m_modelStiffness.push_back(StiffnessType{stiffness});
        m_effectiveStiffness.push_back(StiffnessType{stiffness});
        m_modelCompliance.push_back(m_modelStiffness.back().inverse());
        m_effectiveCompliance.push_back(m_effectiveStiffness.back().inverse());
        m_loadDirection.push_back(loadDirection.normalized());
    }
    const std::vector<Vector3> &loadDirections() const { return m_loadDirection; }

protected:
    virtual void dEnergy(int i, const SimObject &obj, DScalarFirstD &dE) const;
    virtual void dEnergy(int i, const SimObject &obj, DScalarSecondD &dE) const;
    virtual void dConstraints(int i, const SimObject &obj, DConstraintsFirstD &dC) const;
    virtual void dConstraints(int i, const SimObject &obj, DConstraintsSecondD &dC) const;
    std::vector<Vector3> m_loadDirection;
};
}  // namespace spg