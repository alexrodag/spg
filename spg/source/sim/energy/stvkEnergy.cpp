#include <spg/sim/energy/stvkEnergy.h>
#include <spg/sim/simObject.h>
#include <cmath>

namespace spg
{
namespace
{
auto l_stvkConstraint = [](const StvkEnergy *energy, const int i, const SimObject &obj, auto &dC) {
    const auto &x0p{obj.positions()[energy->stencils()[i][0]]};
    const auto &x1p{obj.positions()[energy->stencils()[i][1]]};
    const auto &x2p{obj.positions()[energy->stencils()[i][2]]};
    const auto &x3p{obj.positions()[energy->stencils()[i][3]]};

    using ADouble = std::decay_t<decltype(dC[0])>;
    const Vector3T<ADouble> x0(ADouble(x0p.x(), 0), ADouble(x0p.y(), 1), ADouble(x0p.z(), 2));
    const Vector3T<ADouble> x1(ADouble(x1p.x(), 3), ADouble(x1p.y(), 4), ADouble(x1p.z(), 5));
    const Vector3T<ADouble> x2(ADouble(x2p.x(), 6), ADouble(x2p.y(), 7), ADouble(x2p.z(), 8));
    const Vector3T<ADouble> x3(ADouble(x3p.x(), 9), ADouble(x3p.y(), 10), ADouble(x3p.z(), 11));
    // Compute deformation gradient
    const Vector3T<ADouble> u = x1 - x0;
    const Vector3T<ADouble> v = x2 - x0;
    const Vector3T<ADouble> w = x3 - x0;
    const Matrix3T<ADouble> deformedMatrix{{u[0], v[0], w[0]}, {u[1], v[1], w[1]}, {u[2], v[2], w[2]}};
    const Matrix3T<ADouble> F = deformedMatrix * energy->inverseReferenceMats()[i];
    // Compute green strain
    const Matrix3T<ADouble> E = 0.5 * (F.transpose() * F - Matrix3::Identity());
    const auto hydrostaticTerm = E.trace();
    // add a 2* to compensate the 0.5* introduced by the constraint-based energy computation
    const auto deviatoricTerm = sqrt(2 * (E.transpose() * E).trace());
    dC = {hydrostaticTerm, deviatoricTerm};
    // Note: This constraint version has a numerical instability when the second term is zero, since the derivative of
    // the sqrt would lead to a division by zero. Interestingly, this limitation would not be present in the full
    // analitical version of the energy using the constraint approach, but automatic differentiation is not capable of
    // doing analitical simplifications. Certainly a significant limitation of autodiff in some cases
    if (dC[1].val == 0) {
        dC[1] = 0;
    }
};

auto l_stvkEnergy = [](const StvkEnergy *energy, const int i, const SimObject &obj, auto &dE) {
    const auto &x0p{obj.positions()[energy->stencils()[i][0]]};
    const auto &x1p{obj.positions()[energy->stencils()[i][1]]};
    const auto &x2p{obj.positions()[energy->stencils()[i][2]]};
    const auto &x3p{obj.positions()[energy->stencils()[i][3]]};

    using ADouble = std::decay_t<decltype(dE)>;
    const Vector3T<ADouble> x0(ADouble(x0p.x(), 0), ADouble(x0p.y(), 1), ADouble(x0p.z(), 2));
    const Vector3T<ADouble> x1(ADouble(x1p.x(), 3), ADouble(x1p.y(), 4), ADouble(x1p.z(), 5));
    const Vector3T<ADouble> x2(ADouble(x2p.x(), 6), ADouble(x2p.y(), 7), ADouble(x2p.z(), 8));
    const Vector3T<ADouble> x3(ADouble(x3p.x(), 9), ADouble(x3p.y(), 10), ADouble(x3p.z(), 11));
    // Compute deformation gradient
    const Vector3T<ADouble> u = x1 - x0;
    const Vector3T<ADouble> v = x2 - x0;
    const Vector3T<ADouble> w = x3 - x0;
    const Matrix3T<ADouble> deformedMatrix{{u[0], v[0], w[0]}, {u[1], v[1], w[1]}, {u[2], v[2], w[2]}};
    const Matrix3T<ADouble> F = deformedMatrix * energy->inverseReferenceMats()[i];
    // Compute green strain
    const Matrix3T<ADouble> E = 0.5 * (F.transpose() * F - Matrix3::Identity());
    // Lambda and mu come pre-multiplied by the volume
    const auto lambda = energy->effectiveStiffness()[i](0, 0);
    const auto mu = energy->effectiveStiffness()[i](1, 1);
    const auto hydrostaticTerm = E.trace();
    dE = 0.5 * lambda * hydrostaticTerm * hydrostaticTerm + mu * (E.transpose() * E).trace();
};
}  // namespace

void StvkEnergy::addStencil(const std::array<int, s_stencilSize> &stencil, const Real young, const Real poisson)
{
    m_stencils.push_back(stencil);
    const Real mu = young / (2 * (1 + poisson));
    const Real lambda = -mu / (2 * poisson - 1);

    const Matrix2 C{{lambda, 0}, {0, mu}};
    m_modelStiffness.emplace_back(C);
    m_effectiveStiffness.emplace_back(C);
    m_modelCompliance.emplace_back(C.inverse());
    m_effectiveCompliance.emplace_back(C.inverse());
}

void StvkEnergy::preparePrecomputations(const SimObject &obj)
{
    const int nstencils{static_cast<int>(m_stencils.size())};
    m_inverseReferenceMat.resize(nstencils);
    for (int i = 0; i < nstencils; ++i) {
        const auto &x0{obj.positions0()[m_stencils[i][0]]};
        const auto &x1{obj.positions0()[m_stencils[i][1]]};
        const auto &x2{obj.positions0()[m_stencils[i][2]]};
        const auto &x3{obj.positions0()[m_stencils[i][3]]};
        // compute rest material matrix
        const Vector3 u = x1 - x0;
        const Vector3 v = x2 - x0;
        const Vector3 w = x3 - x0;

        m_inverseReferenceMat[i].col(0) = u;
        m_inverseReferenceMat[i].col(1) = v;
        m_inverseReferenceMat[i].col(2) = w;
        const auto volume = 1. / 6. * m_inverseReferenceMat[i].determinant();
        if (volume <= 0) {
            throw std::runtime_error("SNH tetrahedron rest volume is not positive");
        }
        m_inverseReferenceMat[i] = m_inverseReferenceMat[i].inverse().eval();
        m_effectiveCompliance[i] = m_modelCompliance[i] / volume;
        m_effectiveStiffness[i] = m_modelStiffness[i] * volume;
    }
    StencilBlockEnergy<4, 2>::preparePrecomputations(obj);
}

void StvkEnergy::dConstraints(int i, const SimObject &obj, DConstraintsFirstD &dC) const
{
    l_stvkConstraint(this, i, obj, dC);
}

void StvkEnergy::dConstraints(int i, const SimObject &obj, DConstraintsSecondD &dC) const
{
    l_stvkConstraint(this, i, obj, dC);
}

void StvkEnergy::dEnergy(int i, const SimObject &obj, DScalarFirstD &dE) const
{
    l_stvkEnergy(this, i, obj, dE);
}

void StvkEnergy::dEnergy(int i, const SimObject &obj, DScalarSecondD &dE) const
{
    l_stvkEnergy(this, i, obj, dE);
}
}  // namespace spg