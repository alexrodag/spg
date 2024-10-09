#include <spg/sim/energy/springAnchorEnergy.h>
#include <spg/sim/simObject.h>

namespace spg
{
namespace
{
auto l_springConstraint = [](const SpringAnchorEnergy *energy, const int i, const SimObject &obj, auto &dC) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    using ADouble = std::decay_t<decltype(dC[0])>;
    const Vector3T<ADouble> x(ADouble(x0.x(), 0), ADouble(x0.y(), 1), ADouble(x0.z(), 2));
    dC[0] = (x - energy->anchors()[i]).norm();
    // Required due to undefined gradients when the norm is 0 due to division by 0 in the gradient computation. This is
    // a limitation of using differentiation at the constraint level instead of energy level, where this division doesnt
    // appear
    if (dC[0].val == 0) {
        dC[0] = 0;
    }
};

auto l_springEnergy = [](const SpringAnchorEnergy *energy, const int i, const SimObject &obj, auto &dE) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    using ADouble = std::decay_t<decltype(dE)>;
    const Vector3T<ADouble> x(ADouble(x0.x(), 0), ADouble(x0.y(), 1), ADouble(x0.z(), 2));
    dE = 0.5 * energy->modelStiffness()[i][0] * (x - energy->anchors()[i]).squaredNorm();
};
}  // namespace

void SpringAnchorEnergy::addStencil(std::array<int, s_stencilSize> stencil,
                                    const Vector3 &anchor,
                                    const Real stiffness = 0)
{
    m_stencils.push_back(stencil);
    m_modelStiffness.push_back(StiffnessType{stiffness});
    m_effectiveStiffness.push_back(StiffnessType{stiffness});
    m_modelCompliance.push_back(m_modelStiffness.back().inverse());
    m_effectiveCompliance.push_back(m_effectiveStiffness.back().inverse());
    m_anchor.push_back(anchor);
}

void SpringAnchorEnergy::dEnergy(const int i, const SimObject &obj, DScalarFirstD &dC) const
{
    l_springEnergy(this, i, obj, dC);
}

void SpringAnchorEnergy::dEnergy(const int i, const SimObject &obj, DScalarSecondD &dC) const
{
    l_springEnergy(this, i, obj, dC);
}

void SpringAnchorEnergy::dConstraints(const int i, const SimObject &obj, DConstraintsFirstD &dC) const
{
    l_springConstraint(this, i, obj, dC);
}

void SpringAnchorEnergy::dConstraints(const int i, const SimObject &obj, DConstraintsSecondD &dC) const
{
    l_springConstraint(this, i, obj, dC);
}
}  // namespace spg