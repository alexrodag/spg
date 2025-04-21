#include <spg/sim/energy/rigidBodyOrientationAlignmentEnergy.h>

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

auto l_orientationAlignmentConstraint =
    [](const RBOrientationAlignmentEnergy *energy, const int i, const RigidBodyGroup &rbGroup, auto &dC) {
        using RealT = std::decay_t<decltype(dC[0])>;
        const Vector3T<RealT> theta0(RealT(0, 3), RealT(0, 4), RealT(0, 5));
        const Vector3T<RealT> theta1(RealT(0, 9), RealT(0, 10), RealT(0, 11));
        const Matrix3T<RealT> rotMatrix0 =
            (Matrix3::Identity() + skew(theta0)) * rbGroup.rotationMatrices()[energy->stencils()[i][0]];
        const Matrix3T<RealT> rotMatrix1 =
            (Matrix3::Identity() + skew(theta1)) * rbGroup.rotationMatrices()[energy->stencils()[i][1]];
        const Matrix3T<RealT> deltaRotMatrix = rotMatrix0 * rotMatrix1.transpose();
        dC[0] = rotMatrix0.col(0).cross(rotMatrix1.col(0)).squaredNorm() +
                rotMatrix0.col(1).cross(rotMatrix1.col(1)).squaredNorm() +
                rotMatrix0.col(2).cross(rotMatrix1.col(2)).squaredNorm();
    };
}  // namespace

void RBOrientationAlignmentEnergy::addStencil(std::array<int, s_stencilSize> stencil, const Real stiffness = 0)
{
    m_stencils.push_back(stencil);
    m_modelStiffness.push_back(StiffnessMat{stiffness});
    m_effectiveStiffness.push_back(StiffnessMat{stiffness});
    m_modelCompliance.push_back(m_modelStiffness.back().inverse());
    m_effectiveCompliance.push_back(m_effectiveStiffness.back().inverse());
}

void RBOrientationAlignmentEnergy::dConstraints(const int i, const RigidBodyGroup &rbGroup, ConstraintsAD1 &dC) const
{
    l_orientationAlignmentConstraint(this, i, rbGroup, dC);
}

void RBOrientationAlignmentEnergy::dConstraints(const int i, const RigidBodyGroup &rbGroup, ConstraintsAD2 &dC) const
{
    l_orientationAlignmentConstraint(this, i, rbGroup, dC);
}
}  // namespace spg