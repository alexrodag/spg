#include <spg/sim/energy/springRigidBodyEnergy.h>
#include <spg/sim/simObject/particleGroup.h>

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

auto l_springRBConstraint = [](const SpringRBEnergy *energy, const int i, const RigidBodyGroup &rbGroup, auto &dC) {
    const auto &x0{rbGroup.positions()[energy->stencils()[i][0]]};
    const auto &x1{rbGroup.positions()[energy->stencils()[i][1]]};
    using RealT = std::decay_t<decltype(dC[0])>;
    const Vector3T<RealT> com0(RealT(x0.x(), 0), RealT(x0.y(), 1), RealT(x0.z(), 2));
    const Vector3T<RealT> theta0(RealT(0, 3), RealT(0, 4), RealT(0, 5));
    const Vector3T<RealT> com1(RealT(x1.x(), 6), RealT(x1.y(), 7), RealT(x1.z(), 8));
    const Vector3T<RealT> theta1(RealT(0, 9), RealT(0, 10), RealT(0, 11));
    const Vector3T<RealT> globalRBPoint0 = com0 + (Matrix3::Identity() + skew(theta0)) *
                                                      rbGroup.rotationMatrices()[energy->stencils()[i][0]] *
                                                      energy->localRBPoints()[i][0];
    const Vector3T<RealT> globalRBPoint1 = com1 + (Matrix3::Identity() + skew(theta1)) *
                                                      rbGroup.rotationMatrices()[energy->stencils()[i][1]] *
                                                      energy->localRBPoints()[i][1];
    const Real L0 = energy->restLengths()[i];
    dC[0] = (globalRBPoint0 - globalRBPoint1).norm() - L0;
};
}  // namespace

void SpringRBEnergy::addStencil(std::array<int, s_stencilSize> stencil,
                                const std::array<Vector3, 2> &localRigidBodyPoints,
                                const Real restLength,
                                const Real stiffness = 0)
{
    m_stencils.push_back(stencil);
    if (restLength <= 0) {
        throw std::runtime_error("Zero/negative rest length invalid");
    }
    m_restLength.push_back(restLength);
    m_localRBPoints.push_back(localRigidBodyPoints);
    m_modelStiffness.push_back(StiffnessMat{stiffness});
    m_effectiveStiffness.push_back(StiffnessMat{stiffness});
    m_modelCompliance.push_back(m_modelStiffness.back().inverse());
    m_effectiveCompliance.push_back(m_effectiveStiffness.back().inverse());
}

void SpringRBEnergy::dConstraints(const int i, const RigidBodyGroup &rbGroup, ConstraintsAD1 &dC) const
{
    l_springRBConstraint(this, i, rbGroup, dC);
}

void SpringRBEnergy::dConstraints(const int i, const RigidBodyGroup &rbGroup, ConstraintsAD2 &dC) const
{
    l_springRBConstraint(this, i, rbGroup, dC);
}
}  // namespace spg