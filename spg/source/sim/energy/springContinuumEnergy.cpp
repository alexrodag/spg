#include <spg/sim/energy/springContinuumEnergy.h>
#include <spg/sim/simObject.h>

namespace spg
{
namespace
{
auto l_springConstraint = [](const SpringContinuumEnergy *energy, const int i, const SimObject &obj, auto &dC) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    const auto &x1{obj.positions()[energy->stencils()[i][1]]};
    using ADouble = std::decay_t<decltype(dC[0])>;
    const Vector3T<ADouble> x(ADouble(x0.x(), 0), ADouble(x0.y(), 1), ADouble(x0.z(), 2));
    const Vector3T<ADouble> y(ADouble(x1.x(), 3), ADouble(x1.y(), 4), ADouble(x1.z(), 5));
    const Real L0 = energy->restLengths()[i];
    dC[0] = (x - y).norm() / L0 - 1.;
};
}

void SpringContinuumEnergy::addStencil(const std::array<int, s_stencilSize> &stencil,
                                       const Real restLength,
                                       const Real stiffness)
{
    m_stencils.push_back(stencil);
    if (restLength <= 0) {
        throw std::runtime_error("Zero/negative rest length invalid");
    }
    m_restLength.push_back(restLength);
    m_modelStiffness.push_back(StiffnessType{stiffness});
    // Add integration along element for the effective stiffness
    m_effectiveStiffness.push_back(StiffnessType{stiffness * restLength});
    m_modelCompliance.push_back(m_modelStiffness.back().inverse());
    m_effectiveCompliance.push_back(m_effectiveStiffness.back().inverse());
}

void SpringContinuumEnergy::dConstraints(const int i, const SimObject &obj, DConstraintsFirstD &dC) const
{
    l_springConstraint(this, i, obj, dC);
}

void SpringContinuumEnergy::dConstraints(const int i, const SimObject &obj, DConstraintsSecondD &dC) const
{
    l_springConstraint(this, i, obj, dC);
}
}  // namespace spg