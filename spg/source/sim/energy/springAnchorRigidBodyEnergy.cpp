#include <spg/sim/energy/springAnchorRigidBodyEnergy.h>
#include <spg/sim/simObject.h>

namespace spg
{
namespace
{
template <typename T>
spg::Matrix3T<T> skew(const spg::Vector3T<T> &v)
{
    spg::Matrix3T<T> vSkew;
    vSkew << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
    return vSkew;
}

auto l_springRBConstraint = [](const SpringAnchorRBEnergy *energy, const int i, const RigidBodyGroup &obj, auto &dC) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    using RealT = std::decay_t<decltype(dC[0])>;
    const Vector3T<RealT> x(RealT(x0.x(), 0), RealT(x0.y(), 1), RealT(x0.z(), 2));
    const Vector3T<RealT> theta(RealT(0, 3), RealT(0, 4), RealT(0, 5));
    const Vector3T<RealT> globalRBPoint = x + (Matrix3::Identity() + skew(theta)) *
                                                  obj.rotationMatrices()[energy->stencils()[i][0]] *
                                                  energy->localRBPoints()[i];
    dC[0] = (globalRBPoint - energy->anchors()[i]).norm();
    // Required due to undefined gradients when the norm is 0 due to division by 0 in the gradient computation. This
    // is a limitation of using differentiation at the constraint level instead of energy level, where this division
    // doesnt appear
    if (dC[0].val == 0) {
        dC[0] = 0;
    }
};

auto l_springRBEnergy = [](const SpringAnchorRBEnergy *energy, const int i, const RigidBodyGroup &obj, auto &dE) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    using RealT = std::decay_t<decltype(dE)>;
    const Vector3T<RealT> x(RealT(x0.x(), 0), RealT(x0.y(), 1), RealT(x0.z(), 2));
    const Vector3T<RealT> theta(RealT(0, 3), RealT(0, 4), RealT(0, 5));
    const Vector3T<RealT> globalRBPoint = x + (Matrix3::Identity() + skew(theta)) *
                                                  obj.rotationMatrices()[energy->stencils()[i][0]] *
                                                  energy->localRBPoints()[i];

    dE = 0.5 * energy->modelStiffness()[i][0] * (globalRBPoint - energy->anchors()[i]).squaredNorm();
};

struct RBPointMapping {
    Vector3 m_localPoint;
    Vector3 mapPoint(const Vector3 &x, const Matrix3 &R) { return x + R * m_localPoint; }
    Matrix<3, 6> jacobian(const Matrix3 &R)
    {
        Matrix<3, 6> j;
        j.block<3, 3>(0, 0) = Matrix3::Identity();
        j.block<3, 3>(0, 3) = skew((R * -m_localPoint).eval());
        return j;
    }
};
auto l_springRBConnectorEnergy =
    [](const SpringAnchorRBEnergy *energy, const int i, const RigidBodyGroup &obj, auto &dE) {
        const auto &x{obj.positions()[energy->stencils()[i][0]]};
        const Vector3 globalRBPoint = x + obj.rotationMatrices()[energy->stencils()[i][0]] * energy->localRBPoints()[i];
        using RealT = std::decay_t<decltype(dE)>;
        const Vector3T<RealT> p(RealT(globalRBPoint.x(), 0), RealT(globalRBPoint.y(), 1), RealT(globalRBPoint.z(), 2));
        dE = 0.5 * energy->modelStiffness()[i][0] * (p - energy->anchors()[i]).squaredNorm();
    };
}  // namespace

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

Real SpringAnchorRBEnergy::energy(int i, const RigidBodyGroup &obj) const
{
    TinyAD::Scalar<3, Real, false> e;
    l_springRBConnectorEnergy(this, i, obj, e);
    return e.val;
}
SpringAnchorRBEnergy::EnergyGrad SpringAnchorRBEnergy::energyGradient(int i, const RigidBodyGroup &obj) const
{
    TinyAD::Scalar<3, Real, false> e;
    l_springRBConnectorEnergy(this, i, obj, e);
    RBPointMapping mapping;
    mapping.m_localPoint = this->localRBPoints()[i];
    const auto J = mapping.jacobian(obj.rotationMatrices()[i]);
    return J.transpose() * e.grad;
}
SpringAnchorRBEnergy::EnergyHess SpringAnchorRBEnergy::energyHessian(int i, const RigidBodyGroup &obj) const
{
    TinyAD::Scalar<3, Real> e;
    l_springRBConnectorEnergy(this, i, obj, e);
    RBPointMapping mapping;
    mapping.m_localPoint = this->localRBPoints()[i];
    const auto J = mapping.jacobian(obj.rotationMatrices()[i]);
    return J.transpose() * e.Hess * J;
}
}  // namespace spg