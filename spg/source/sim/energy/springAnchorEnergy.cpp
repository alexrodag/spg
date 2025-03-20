#include <spg/sim/energy/springAnchorEnergy.h>
#include <spg/sim/simObject/particleGroup.h>

namespace spg
{
namespace
{
auto l_springConstraint = [](const SpringAnchorEnergy *energy, const int i, const ParticleGroup &pGroup, auto &dC) {
    const auto &x0{pGroup.positions()[energy->stencils()[i][0]]};
    using RealT = std::decay_t<decltype(dC[0])>;
    const Vector3T<RealT> x(RealT(x0.x(), 0), RealT(x0.y(), 1), RealT(x0.z(), 2));
    dC[0] = (x - energy->anchors()[i]).norm();
    // Required due to undefined gradients when the norm is 0 due to division by 0 in the gradient computation. This is
    // a limitation of using differentiation at the constraint level instead of energy level, where this division doesnt
    // appear
    if (dC[0].val == 0) {
        dC[0] = 0;
    }
};

auto l_springEnergy = [](const SpringAnchorEnergy *energy, const int i, const ParticleGroup &pGroup, auto &dE) {
    const auto &x0{pGroup.positions()[energy->stencils()[i][0]]};
    using RealT = std::decay_t<decltype(dE)>;
    const Vector3T<RealT> x(RealT(x0.x(), 0), RealT(x0.y(), 1), RealT(x0.z(), 2));
    dE = 0.5 * energy->modelStiffness()[i][0] * (x - energy->anchors()[i]).squaredNorm();
};
}  // namespace

void SpringAnchorEnergy::addStencil(std::array<int, s_stencilSize> stencil,
                                    const Vector3 &anchor,
                                    const Real stiffness = 0)
{
    m_stencils.push_back(stencil);
    m_modelStiffness.push_back(StiffnessMat{stiffness});
    m_effectiveStiffness.push_back(StiffnessMat{stiffness});
    m_modelCompliance.push_back(m_modelStiffness.back().inverse());
    m_effectiveCompliance.push_back(m_effectiveStiffness.back().inverse());
    m_anchor.push_back(anchor);
}

void SpringAnchorEnergy::dEnergy(const int i, const ParticleGroup &pGroup, RealAD1 &dC) const
{
    l_springEnergy(this, i, pGroup, dC);
}

void SpringAnchorEnergy::dEnergy(const int i, const ParticleGroup &pGroup, RealAD2 &dC) const
{
    l_springEnergy(this, i, pGroup, dC);
}

void SpringAnchorEnergy::dConstraints(const int i, const ParticleGroup &pGroup, ConstraintsAD1 &dC) const
{
    l_springConstraint(this, i, pGroup, dC);
}

void SpringAnchorEnergy::dConstraints(const int i, const ParticleGroup &pGroup, ConstraintsAD2 &dC) const
{
    l_springConstraint(this, i, pGroup, dC);
}

}  // namespace spg