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

    using RealT = std::decay_t<decltype(dC[0])>;
    const Vector3T<RealT> x0(RealT(x0p.x(), 0), RealT(x0p.y(), 1), RealT(x0p.z(), 2));
    const Vector3T<RealT> x1(RealT(x1p.x(), 3), RealT(x1p.y(), 4), RealT(x1p.z(), 5));
    const Vector3T<RealT> x2(RealT(x2p.x(), 6), RealT(x2p.y(), 7), RealT(x2p.z(), 8));
    const Vector3T<RealT> x3(RealT(x3p.x(), 9), RealT(x3p.y(), 10), RealT(x3p.z(), 11));
    // Compute deformation gradient
    const Vector3T<RealT> u = x1 - x0;
    const Vector3T<RealT> v = x2 - x0;
    const Vector3T<RealT> w = x3 - x0;
    const Matrix3T<RealT> deformedMatrix{{u[0], v[0], w[0]}, {u[1], v[1], w[1]}, {u[2], v[2], w[2]}};
    const Matrix3T<RealT> F = deformedMatrix * energy->inverseReferenceMats()[i];
    // Compute green strain
    const Matrix3T<RealT> E = 0.5 * (F.transpose() * F - Matrix3::Identity());
    const auto hydrostaticTerm = E.trace();
    // add a 2* to compensate the 0.5* introduced by the constraint-based energy computation
    const auto deviatoricTerm = sqrt(2 * (E * E).trace());
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

    using RealT = std::decay_t<decltype(dE)>;
    const Vector3T<RealT> x0(RealT(x0p.x(), 0), RealT(x0p.y(), 1), RealT(x0p.z(), 2));
    const Vector3T<RealT> x1(RealT(x1p.x(), 3), RealT(x1p.y(), 4), RealT(x1p.z(), 5));
    const Vector3T<RealT> x2(RealT(x2p.x(), 6), RealT(x2p.y(), 7), RealT(x2p.z(), 8));
    const Vector3T<RealT> x3(RealT(x3p.x(), 9), RealT(x3p.y(), 10), RealT(x3p.z(), 11));
    // Compute deformation gradient
    const Vector3T<RealT> u = x1 - x0;
    const Vector3T<RealT> v = x2 - x0;
    const Vector3T<RealT> w = x3 - x0;
    const Matrix3T<RealT> deformedMatrix{{u[0], v[0], w[0]}, {u[1], v[1], w[1]}, {u[2], v[2], w[2]}};
    const Matrix3T<RealT> F = deformedMatrix * energy->inverseReferenceMats()[i];
    // Compute green strain
    const Matrix3T<RealT> E = 0.5 * (F.transpose() * F - Matrix3::Identity());
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
    // Compute lamé parameters: https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
    const Real mu = young / (2 * (1 + poisson));
    // Add epsilon to lambda just in case it becomes 0, to avoid division by 0 when computing the compliance
    constexpr Real epsilon = 1e-20;
    const Real lambda = young * poisson / ((1 + poisson) * (1 - 2 * poisson)) + epsilon;

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

void StvkEnergy::dConstraints(int i, const SimObject &obj, ConstraintsAD1 &dC) const
{
    l_stvkConstraint(this, i, obj, dC);
}

void StvkEnergy::dConstraints(int i, const SimObject &obj, ConstraintsAD2 &dC) const
{
    l_stvkConstraint(this, i, obj, dC);
}

void StvkEnergy::dEnergy(int i, const SimObject &obj, RealAD1 &dE) const
{
    l_stvkEnergy(this, i, obj, dE);
}

void StvkEnergy::dEnergy(int i, const SimObject &obj, RealAD2 &dE) const
{
    l_stvkEnergy(this, i, obj, dE);
}
}  // namespace spg