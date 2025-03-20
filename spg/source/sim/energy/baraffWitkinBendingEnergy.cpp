#include <spg/sim/energy/baraffWitkinBendingEnergy.h>
#include <spg/sim/simObject/particleGroup.h>
#include <cmath>

namespace spg
{
namespace
{
auto l_bendingConstraint =
    [](const BaraffWitkinBendingEnergy *energy, const int i, const ParticleGroup &obj, auto &dC) {
        const auto &x0p{obj.positions()[energy->stencils()[i][0]]};
        const auto &x1p{obj.positions()[energy->stencils()[i][1]]};
        const auto &x2p{obj.positions()[energy->stencils()[i][2]]};
        const auto &x3p{obj.positions()[energy->stencils()[i][3]]};
        using RealT = std::decay_t<decltype(dC[0])>;
        const Vector3T<RealT> x0(RealT(x0p.x(), 0), RealT(x0p.y(), 1), RealT(x0p.z(), 2));
        const Vector3T<RealT> x1(RealT(x1p.x(), 3), RealT(x1p.y(), 4), RealT(x1p.z(), 5));
        const Vector3T<RealT> x2(RealT(x2p.x(), 6), RealT(x2p.y(), 7), RealT(x2p.z(), 8));
        const Vector3T<RealT> x3(RealT(x3p.x(), 9), RealT(x3p.y(), 10), RealT(x3p.z(), 11));

        const Vector3T<RealT> e0 = x1 - x0;
        const Vector3T<RealT> e3 = x2 - x1;
        const Vector3T<RealT> e4 = x3 - x1;

        // use non-differentiable norms to avoid derivatives w.r.t. them, as suggested in the paper
        const Real e0NormInv = 1 / e0.norm().val;

        const Vector3T<RealT> n1 = e0.cross(e3);
        const Real n1NormInv = 1 / n1.norm().val;
        const Vector3T<RealT> n2 = -e0.cross(e4);
        const Real n2NormInv = 1 / n2.norm().val;

        const auto theta =
            atan2(((n1 * n1NormInv).cross(n2 * n2NormInv)).dot(e0 * e0NormInv), (n1 * n1NormInv).dot(n2 * n2NormInv));

        dC[0] = theta - energy->restThetas()[i];
    };
}

void BaraffWitkinBendingEnergy::addStencil(const std::array<int, s_stencilSize> &stencil,
                                           const Real restTheta,
                                           const Real stiffness)
{
    m_stencils.push_back(stencil);
    if (restTheta < -3.12 || restTheta > 3.12) {
        throw std::runtime_error("Invalid restTheta");
    }
    m_restTheta.push_back(restTheta);
    m_modelStiffness.push_back(StiffnessMat{stiffness});
    m_effectiveStiffness.push_back(StiffnessMat{stiffness});
    m_modelCompliance.push_back(m_modelStiffness.back().inverse());
    m_effectiveCompliance.push_back(m_effectiveStiffness.back().inverse());
}

void BaraffWitkinBendingEnergy::preparePrecomputations(const ParticleGroup &obj)
{
    const int nstencils{static_cast<int>(m_stencils.size())};
    for (int i = 0; i < nstencils; ++i) {
        const auto &x0{obj.positions0()[m_stencils[i][0]]};
        const auto &x1{obj.positions0()[m_stencils[i][1]]};
        const auto &x2{obj.positions0()[m_stencils[i][2]]};
        const auto &x3{obj.positions0()[m_stencils[i][3]]};

        const Vector3 e0 = x1 - x0;
        const Vector3 e3 = x2 - x1;
        const Vector3 e4 = x3 - x1;

        // In contrast to normal dihedral bendings, the paper suggest not using area, to better match the concept of
        // curvature instead of angle

        // Compute rest angles based on initial positions
        const Vector3 n1 = e0.cross(e3).normalized();
        const Vector3 n2 = -e0.cross(e4).normalized();
        m_restTheta[i] = atan2((n1.cross(n2)).dot(e0.normalized()), n1.dot(n2));
    }
    StencilBlockEnergy<4>::preparePrecomputations(obj);
}

void BaraffWitkinBendingEnergy::dConstraints(int i, const ParticleGroup &obj, ConstraintsAD1 &dC) const
{
    l_bendingConstraint(this, i, obj, dC);
}
void BaraffWitkinBendingEnergy::dConstraints(int i, const ParticleGroup &obj, ConstraintsAD2 &dC) const
{
    l_bendingConstraint(this, i, obj, dC);
}

}  // namespace spg