#include <spg/sim/energy/stableNeoHookeanEnergy.h>
#include <spg/sim/simObject/particleGroup.h>
#include <cmath>

namespace spg
{
namespace
{
auto l_stableNeoHookeanConstraint =
    [](const StableNeoHookeanEnergy *energy, const int i, const ParticleGroup &obj, auto &dC) {
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
        // Compute constraint terms
        const auto hydrostaticTerm = F.determinant() - energy->hydrostaticCorrectionTerms()[i];
        const auto deviatoricTerm = sqrt((F.transpose() * F).trace());
        dC = {hydrostaticTerm, deviatoricTerm};
    };

auto l_stableNeoHookeanEnergy =
    [](const StableNeoHookeanEnergy *energy, const int i, const ParticleGroup &obj, auto &dE) {
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
        // Compute energy
        // Lambda and mu come pre-multiplied by the volume
        // Note: This implementation could be removed, as the automatic constraint-based energy computation based on the
        // constraint implementation above would produce equivalent forces and hessians. This implementation is mainly
        // for didactic purposes
        const auto lambda = energy->effectiveStiffness()[i](0, 0);
        const auto mu = energy->effectiveStiffness()[i](1, 1);
        const auto hydrostaticTerm = F.determinant() - energy->hydrostaticCorrectionTerms()[i];
        const auto deviatoricTerm = (F.transpose() * F).trace() - 3;
        dE = 0.5 * lambda * hydrostaticTerm * hydrostaticTerm + 0.5 * mu * deviatoricTerm;
    };
}  // namespace

void StableNeoHookeanEnergy::addStencil(const std::array<int, s_stencilSize> &stencil,
                                        const Real young,
                                        const Real poisson)
{
    m_stencils.push_back(stencil);
    // Note: These formulas are derived from the correspondences mu = mu_lamé and lambda = mu_lamé + lambda_lamé as
    // explained in sec 3.4 of Smith 2018
    const Real mu = young / (2 * (1 + poisson));
    const Real lambda = mu / (1 - 2 * poisson);

    const Matrix2 C{{lambda, 0}, {0, mu}};
    m_modelStiffness.emplace_back(C);
    m_effectiveStiffness.emplace_back(C);
    m_modelCompliance.emplace_back(C.inverse());
    m_effectiveCompliance.emplace_back(C.inverse());
    m_hydrostaticCorrectionTerm.emplace_back(1 + mu / lambda);
}

void StableNeoHookeanEnergy::preparePrecomputations(const ParticleGroup &obj)
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

void StableNeoHookeanEnergy::dConstraints(int i, const ParticleGroup &obj, ConstraintsAD1 &dC) const
{
    l_stableNeoHookeanConstraint(this, i, obj, dC);
}

void StableNeoHookeanEnergy::dConstraints(int i, const ParticleGroup &obj, ConstraintsAD2 &dC) const
{
    l_stableNeoHookeanConstraint(this, i, obj, dC);
}

void StableNeoHookeanEnergy::dEnergy(int i, const ParticleGroup &obj, RealAD1 &dE) const
{
    l_stableNeoHookeanEnergy(this, i, obj, dE);
}

void StableNeoHookeanEnergy::dEnergy(int i, const ParticleGroup &obj, RealAD2 &dE) const
{
    l_stableNeoHookeanEnergy(this, i, obj, dE);
}
}  // namespace spg