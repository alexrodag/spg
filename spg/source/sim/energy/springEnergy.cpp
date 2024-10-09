#include <spg/sim/energy/springEnergy.h>
#include <spg/sim/simObject.h>

namespace spg
{
namespace
{
constexpr Real eps{static_cast<Real>(1.0e-6)};  // TODO decide if this should be moved somewhere else

inline Vector3 getDistanceCorrection(const Vector3 &p0,
                                     const Vector3 &p1,
                                     const Real restLength,
                                     const Real alphaTilde,
                                     const Real wSum)
{
    Vector3 n{p1 - p0};
    const Real d{n.norm()};
    if (d < eps) {
        return {0, 0, 0};
    }
    n = n / d;

    return n * (d - restLength) / (wSum + alphaTilde);
}

auto l_springConstraint = [](const SpringEnergy *energy, const int i, const SimObject &obj, auto &dC) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    const auto &x1{obj.positions()[energy->stencils()[i][1]]};
    using ADouble = std::decay_t<decltype(dC[0])>;
    const Vector3T<ADouble> x(ADouble(x0.x(), 0), ADouble(x0.y(), 1), ADouble(x0.z(), 2));
    const Vector3T<ADouble> y(ADouble(x1.x(), 3), ADouble(x1.y(), 4), ADouble(x1.z(), 5));
    const Real L0 = energy->restLengths()[i];
    dC[0] = (x - y).norm() - L0;
};

auto l_springEnergy = [](const SpringEnergy *energy, const int i, const SimObject &obj, auto &dE) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    const auto &x1{obj.positions()[energy->stencils()[i][1]]};
    using ADouble = std::decay_t<decltype(dE)>;
    const Vector3T<ADouble> x(ADouble(x0.x(), 0), ADouble(x0.y(), 1), ADouble(x0.z(), 2));
    const Vector3T<ADouble> y(ADouble(x1.x(), 3), ADouble(x1.y(), 4), ADouble(x1.z(), 5));
    const Real k = energy->effectiveStiffness()[i][0];
    const Real L0 = energy->restLengths()[i];
    const auto diff = (x - y).norm() - L0;
    dE = 0.5 * k * diff * diff;
};
}  // namespace

void SpringEnergy::addStencil(const std::array<int, s_stencilSize> &stencil,
                              const Real restLength,
                              const Real stiffness)
{
    m_stencils.push_back(stencil);
    if (restLength <= 0) {
        throw std::runtime_error("Zero/negative rest length invalid");
    }
    m_restLength.push_back(restLength);
    m_modelStiffness.push_back(StiffnessType{stiffness});
    m_effectiveStiffness.push_back(StiffnessType{stiffness});
    m_modelCompliance.push_back(m_modelStiffness.back().inverse());
    m_effectiveCompliance.push_back(m_effectiveStiffness.back().inverse());
}

void SpringEnergy::projectPosition(const int i, SimObject &obj, const Real dt) const
{
    auto &x0{obj.positions()[m_stencils[i][0]]};
    auto &x1{obj.positions()[m_stencils[i][1]]};
    const auto invMass0{obj.invMasses()[m_stencils[i][0]]};
    const auto invMass1{obj.invMasses()[m_stencils[i][1]]};
    const Real wSum{invMass0 + invMass1};
    if (wSum < eps) {
        return;
    }
    const auto alphaTilde = m_effectiveCompliance[i][0] / (dt * dt);
    Vector3 n{x1 - x0};
    const Real d{n.norm()};
    if (d < eps) {
        return;
    }
    n = n / d;
    const auto corr = n * (d - m_restLength[i]) / (wSum + alphaTilde);
    x0 += invMass0 * corr;
    x1 -= invMass1 * corr;
}

void SpringEnergy::dEnergy(const int i, const SimObject &obj, DScalarFirstD &dC) const
{
    l_springEnergy(this, i, obj, dC);
}

void SpringEnergy::dEnergy(const int i, const SimObject &obj, DScalarSecondD &dC) const
{
    l_springEnergy(this, i, obj, dC);
}

void SpringEnergy::dConstraints(const int i, const SimObject &obj, DConstraintsFirstD &dC) const
{
    l_springConstraint(this, i, obj, dC);
}

void SpringEnergy::dConstraints(const int i, const SimObject &obj, DConstraintsSecondD &dC) const
{
    l_springConstraint(this, i, obj, dC);
}

Real SpringEnergy::energy(const int i, const SimObject &obj) const
{
    auto &x0{obj.positions()[m_stencils[i][0]]};
    auto &x1{obj.positions()[m_stencils[i][1]]};

    const Real k = m_effectiveStiffness[i][0];
    const Real L0 = m_restLength[i];
    const auto LTerm = (x0 - x1).norm() - L0;
    return 0.5 * k * LTerm * LTerm;
}

SpringEnergy::EnergyGradType SpringEnergy::energyGradient(const int i, const SimObject &obj) const
{
    const auto &x0{obj.positions()[m_stencils[i][0]]};
    const auto &x1{obj.positions()[m_stencils[i][1]]};

    const Real k = m_effectiveStiffness[i][0];
    const Real L0 = m_restLength[i];
    SpringEnergy::EnergyGradType grad;
    const Vector3 s = x0 - x1;
    const Real sNorm = s.norm();
    const Vector3 gradx0 = k * (sNorm - L0) * s / sNorm;
    grad.segment<3>(0) = gradx0;
    grad.segment<3>(3) = -gradx0;
    return grad;
}

SpringEnergy::EnergyHessType SpringEnergy::energyHessian(const int i, const SimObject &obj) const
{
    const auto &x0{obj.positions()[m_stencils[i][0]]};
    const auto &x1{obj.positions()[m_stencils[i][1]]};

    const Real k = m_effectiveStiffness[i][0];
    const Real L0 = m_restLength[i];
    SpringEnergy::EnergyHessType hess;
    const Vector3 s = x0 - x1;
    const Real sNorm = s.norm();
    const Vector3 sNormalized = s / sNorm;
    const Matrix3 sst = sNormalized * sNormalized.transpose();
    hess.block<3, 3>(0, 0) = k * ((1 - L0 / sNorm) * (Matrix3::Identity() - sst) + sst);
    hess.block<3, 3>(3, 3) = hess.block<3, 3>(0, 0);
    hess.block<3, 3>(0, 3) = -hess.block<3, 3>(0, 0);
    hess.block<3, 3>(3, 0) = hess.block<3, 3>(0, 3);
    return hess;
}
}  // namespace spg