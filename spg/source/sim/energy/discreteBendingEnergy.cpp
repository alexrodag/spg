#include <spg/sim/energy/discreteBendingEnergy.h>
#include <spg/sim/simObject.h>
#include <cmath>

namespace spg
{
namespace
{
auto l_bendingConstraint = [](const DiscreteBendingEnergy *energy, const int i, const SimObject &obj, auto &dC) {
    const auto &x0p{obj.positions()[energy->stencils()[i][0]]};
    const auto &x1p{obj.positions()[energy->stencils()[i][1]]};
    const auto &x2p{obj.positions()[energy->stencils()[i][2]]};
    const auto &x3p{obj.positions()[energy->stencils()[i][3]]};
    using ADouble = std::decay_t<decltype(dC[0])>;
    const Vector3T<ADouble> x0(ADouble(x0p.x(), 0), ADouble(x0p.y(), 1), ADouble(x0p.z(), 2));
    const Vector3T<ADouble> x1(ADouble(x1p.x(), 3), ADouble(x1p.y(), 4), ADouble(x1p.z(), 5));
    const Vector3T<ADouble> x2(ADouble(x2p.x(), 6), ADouble(x2p.y(), 7), ADouble(x2p.z(), 8));
    const Vector3T<ADouble> x3(ADouble(x3p.x(), 9), ADouble(x3p.y(), 10), ADouble(x3p.z(), 11));

    const Vector3T<ADouble> e0 = x1 - x0;
    const Vector3T<ADouble> e3 = x2 - x1;
    const Vector3T<ADouble> e4 = x3 - x1;

    const Vector3T<ADouble> n1 = e0.cross(e3).normalized();
    const Vector3T<ADouble> n2 = -e0.cross(e4).normalized();

    const auto theta = atan2((n1.cross(n2)).dot(e0.normalized()), n1.dot(n2));

    dC[0] = theta - energy->restThetas()[i];
};
}

void DiscreteBendingEnergy::addStencil(const std::array<int, s_stencilSize> &stencil,
                                       const Real restTheta,
                                       const Real stiffness)
{
    m_stencils.push_back(stencil);
    if (restTheta < -3.12 || restTheta > 3.12) {
        throw std::runtime_error("Invalid restTheta");
    }
    m_restTheta.push_back(restTheta);
    m_modelStiffness.push_back(StiffnessType{stiffness});
    m_effectiveStiffness.push_back(StiffnessType{stiffness});
    m_modelCompliance.push_back(m_modelStiffness.back().inverse());
    m_effectiveCompliance.push_back(m_effectiveStiffness.back().inverse());
}

void DiscreteBendingEnergy::preparePrecomputations(const SimObject &obj)
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
        // 3 * ||e_0||^2 / A, as per the paper. Adding a factor of 2 to compensate the 1/2 introduced in the energy
        // computation.
        const auto constantFactor = 2. * 3. * e0.squaredNorm() / (0.5 * (e0.cross(e3).norm() + e0.cross(e4).norm()));
        m_effectiveCompliance[i] = m_modelCompliance[i] / constantFactor;
        m_effectiveStiffness[i] = m_effectiveCompliance[i].inverse();

        // Compute rest angles based on initial positions
        const Vector3 n1 = e0.cross(e3).normalized();
        const Vector3 n2 = -e0.cross(e4).normalized();
        m_restTheta[i] = atan2((n1.cross(n2)).dot(e0.normalized()), n1.dot(n2));
    }
    StencilBlockEnergy<4>::preparePrecomputations(obj);
}

void DiscreteBendingEnergy::dConstraints(int i, const SimObject &obj, DConstraintsFirstD &dC) const
{
    l_bendingConstraint(this, i, obj, dC);
}
void DiscreteBendingEnergy::dConstraints(int i, const SimObject &obj, DConstraintsSecondD &dC) const
{
    l_bendingConstraint(this, i, obj, dC);
}

}  // namespace spg