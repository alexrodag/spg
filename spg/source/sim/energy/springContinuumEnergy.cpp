#include <spg/sim/energy/springContinuumEnergy.h>
#include <spg/sim/simObject/particleGroup.h>

namespace spg
{
namespace
{
auto l_springConstraint = [](const SpringContinuumEnergy *energy, const int i, const ParticleGroup &obj, auto &dC) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    const auto &x1{obj.positions()[energy->stencils()[i][1]]};
    using RealT = std::decay_t<decltype(dC[0])>;
    const Vector3T<RealT> x(RealT(x0.x(), 0), RealT(x0.y(), 1), RealT(x0.z(), 2));
    const Vector3T<RealT> y(RealT(x1.x(), 3), RealT(x1.y(), 4), RealT(x1.z(), 5));
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
    m_modelStiffness.push_back(StiffnessMat{stiffness});
    // Add integration along element for the effective stiffness
    m_effectiveStiffness.push_back(StiffnessMat{stiffness * restLength});
    m_modelCompliance.push_back(m_modelStiffness.back().inverse());
    m_effectiveCompliance.push_back(m_effectiveStiffness.back().inverse());
}

void SpringContinuumEnergy::dConstraints(const int i, const ParticleGroup &obj, ConstraintsAD1 &dC) const
{
    l_springConstraint(this, i, obj, dC);
}

void SpringContinuumEnergy::dConstraints(const int i, const ParticleGroup &obj, ConstraintsAD2 &dC) const
{
    l_springConstraint(this, i, obj, dC);
}
}  // namespace spg