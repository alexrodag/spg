#include <spg/sim/energy/quadraticBendingEnergy.h>
#include <spg/sim/simObject/particleGroup.h>
#include <cmath>

namespace spg
{
namespace
{
auto l_bendingConstraint = [](const QuadraticBendingEnergy *energy, const int i, const ParticleGroup &obj, auto &dC) {
    const auto &x0p{obj.positions()[energy->stencils()[i][0]]};
    const auto &x1p{obj.positions()[energy->stencils()[i][1]]};
    const auto &x2p{obj.positions()[energy->stencils()[i][2]]};
    const auto &x3p{obj.positions()[energy->stencils()[i][3]]};
    using RealT = std::decay_t<decltype(dC[0])>;
    const VectorT<RealT, 12> x(RealT(x0p.x(), 0),
                               RealT(x0p.y(), 1),
                               RealT(x0p.z(), 2),
                               RealT(x1p.x(), 3),
                               RealT(x1p.y(), 4),
                               RealT(x1p.z(), 5),
                               RealT(x2p.x(), 6),
                               RealT(x2p.y(), 7),
                               RealT(x2p.z(), 8),
                               RealT(x3p.x(), 9),
                               RealT(x3p.y(), 10),
                               RealT(x3p.z(), 11));
    dC = (energy->laplacians()[i] * x).eval();
};
}

void QuadraticBendingEnergy::addStencil(const std::array<int, s_stencilSize> &stencil,
                                        const Real restTheta,
                                        const Real stiffness)
{
    // TODO: Check if its possible to add support for this
    if (restTheta != 0) {
        throw std::runtime_error("Rest theta != 0 not supported for quadratic bending");
    }
    m_stencils.push_back(stencil);
    m_laplacians.push_back(Matrix<3, 12>{});
    m_stiffnessFactor.push_back(stiffness);
    m_modelCompliance.push_back(StiffnessMat{});
    m_effectiveCompliance.push_back(StiffnessMat{});
    m_modelStiffness.push_back(StiffnessMat{});
    m_effectiveStiffness.push_back(StiffnessMat{});
}

void QuadraticBendingEnergy::preparePrecomputations(const ParticleGroup &obj)
{
    const int nstencils{static_cast<int>(m_stencils.size())};
    for (int i = 0; i < nstencils; ++i) {
        const auto &x0{obj.positions0()[m_stencils[i][0]]};
        const auto &x1{obj.positions0()[m_stencils[i][1]]};
        const auto &x2{obj.positions0()[m_stencils[i][2]]};
        const auto &x3{obj.positions0()[m_stencils[i][3]]};

        const Vector3 e0 = x1 - x0;
        const Vector3 e1 = x2 - x0;
        const Vector3 e2 = x3 - x0;
        const Vector3 e3 = x2 - x1;
        const Vector3 e4 = x3 - x1;

        auto l_cotTheta = [](const Vector3 &a, const Vector3 &b) {
            const auto cosTheta = a.dot(b);
            const auto sinTheta = a.cross(b).norm();
            return cosTheta / sinTheta;
        };

        const Real c01 = l_cotTheta(e0, e1);
        const Real c02 = l_cotTheta(e0, e2);
        const Real c03 = l_cotTheta(-e0, e3);
        const Real c04 = l_cotTheta(-e0, e4);

        const Real c0 = c03 + c04;
        const Real c1 = c01 + c02;
        const Real c2 = -c01 - c03;
        const Real c3 = -c02 - c04;
        Matrix<12, 3> K;
        K.block<3, 3>(0, 0) = c0 * Matrix3::Identity();
        K.block<3, 3>(3, 0) = c1 * Matrix3::Identity();
        K.block<3, 3>(6, 0) = c2 * Matrix3::Identity();
        K.block<3, 3>(9, 0) = c3 * Matrix3::Identity();
        const StiffnessMat M_inverse =
            m_stiffnessFactor[i] * StiffnessMat::Identity() * (3 / (0.5 * (e0.cross(e1).norm() + e0.cross(e2).norm())));
        m_laplacians[i] = K.transpose();
        m_effectiveStiffness[i] = M_inverse;
        m_modelStiffness[i] = m_effectiveStiffness[i];
        m_effectiveCompliance[i] = M_inverse.inverse();
        m_modelCompliance[i] = m_effectiveCompliance[i];
    }
    StencilBlockEnergy<4, 3>::preparePrecomputations(obj);
}

void QuadraticBendingEnergy::dConstraints(int i, const ParticleGroup &obj, ConstraintsAD1 &dC) const
{
    l_bendingConstraint(this, i, obj, dC);
}
void QuadraticBendingEnergy::dConstraints(int i, const ParticleGroup &obj, ConstraintsAD2 &dC) const
{
    l_bendingConstraint(this, i, obj, dC);
}

}  // namespace spg