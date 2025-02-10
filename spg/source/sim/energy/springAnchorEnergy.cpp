#include <spg/sim/energy/springAnchorEnergy.h>
#include <spg/sim/simObject.h>

namespace spg
{
namespace
{
auto l_springConstraint = [](const SpringAnchorEnergy *energy, const int i, const SimObject &obj, auto &dC) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
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

auto l_springEnergy = [](const SpringAnchorEnergy *energy, const int i, const SimObject &obj, auto &dE) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    using RealT = std::decay_t<decltype(dE)>;
    const Vector3T<RealT> x(RealT(x0.x(), 0), RealT(x0.y(), 1), RealT(x0.z(), 2));
    dE = 0.5 * energy->modelStiffness()[i][0] * (x - energy->anchors()[i]).squaredNorm();
};

// RB

auto l_springRBConstraint = [](const SpringAnchorRBEnergy *energy, const int i, const RigidBodyGroup &obj, auto &dC) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    using RealT = std::decay_t<decltype(dC[0])>;
    const Vector3T<RealT> x(RealT(x0.x(), 0), RealT(x0.y(), 1), RealT(x0.z(), 2));
    const Vector3T<RealT> theta(RealT(0, 3), RealT(0, 4), RealT(0, 5));

    auto l_skew = [](const spg::Vector3T<RealT> &v) {
        spg::Matrix3T<RealT> vSkew;
        vSkew << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
        return vSkew;
    };
    const Vector3T<RealT> globalRBPoint = x + (Matrix3::Identity() + l_skew(theta)) *
                                                  obj.rotationMatrices()[energy->stencils()[i][0]] *
                                                  energy->localRBPoints()[i];
    dC[0] = (globalRBPoint - energy->anchors()[i]).norm();
    // Required due to undefined gradients when the norm is 0 due to division by 0 in the gradient computation. This is
    // a limitation of using differentiation at the constraint level instead of energy level, where this division doesnt
    // appear
    if (dC[0].val == 0) {
        dC[0] = 0;
    }
};

auto l_springRBEnergy = [](const SpringAnchorRBEnergy *energy, const int i, const RigidBodyGroup &obj, auto &dE) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    using RealT = std::decay_t<decltype(dE)>;
    const Vector3T<RealT> x(RealT(x0.x(), 0), RealT(x0.y(), 1), RealT(x0.z(), 2));
    const Vector3T<RealT> theta(RealT(0, 3), RealT(0, 4), RealT(0, 5));

    auto l_skew = [](const spg::Vector3T<RealT> &v) {
        spg::Matrix3T<RealT> vSkew;
        vSkew << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
        return vSkew;
    };
    const Vector3T<RealT> globalRBPoint = x + (Matrix3::Identity() + l_skew(theta)) *
                                                  obj.rotationMatrices()[energy->stencils()[i][0]] *
                                                  energy->localRBPoints()[i];

    dE = 0.5 * energy->modelStiffness()[i][0] * (globalRBPoint - energy->anchors()[i]).squaredNorm();
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

void SpringAnchorEnergy::dEnergy(const int i, const SimObject &obj, RealAD1 &dC) const
{
    l_springEnergy(this, i, obj, dC);
}

void SpringAnchorEnergy::dEnergy(const int i, const SimObject &obj, RealAD2 &dC) const
{
    l_springEnergy(this, i, obj, dC);
}

void SpringAnchorEnergy::dConstraints(const int i, const SimObject &obj, ConstraintsAD1 &dC) const
{
    l_springConstraint(this, i, obj, dC);
}

void SpringAnchorEnergy::dConstraints(const int i, const SimObject &obj, ConstraintsAD2 &dC) const
{
    l_springConstraint(this, i, obj, dC);
}

// RB

void SpringAnchorRBEnergy::addStencil(std::array<int, s_stencilSize> stencil,
                                      const Vector3 &anchor,
                                      const Vector3 &localRigidBodyPoint,
                                      const Real stiffness = 0)
{
    m_stencils.push_back(stencil);
    m_modelStiffness.push_back(StiffnessMat{stiffness});
    m_effectiveStiffness.push_back(StiffnessMat{stiffness});
    m_modelCompliance.push_back(m_modelStiffness.back().inverse());
    m_effectiveCompliance.push_back(m_effectiveStiffness.back().inverse());
    m_anchor.push_back(anchor);
    m_localRBPoint.push_back(localRigidBodyPoint);
}

void SpringAnchorRBEnergy::dEnergy(const int i, const RigidBodyGroup &obj, RealAD1 &dC) const
{
    l_springRBEnergy(this, i, obj, dC);
}

void SpringAnchorRBEnergy::dEnergy(const int i, const RigidBodyGroup &obj, RealAD2 &dC) const
{
    l_springRBEnergy(this, i, obj, dC);
}

void SpringAnchorRBEnergy::dConstraints(const int i, const RigidBodyGroup &obj, ConstraintsAD1 &dC) const
{
    l_springRBConstraint(this, i, obj, dC);
}

void SpringAnchorRBEnergy::dConstraints(const int i, const RigidBodyGroup &obj, ConstraintsAD2 &dC) const
{
    l_springRBConstraint(this, i, obj, dC);
}
}  // namespace spg