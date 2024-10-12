#include <spg/sim/energy/membraneStvkEnergy.h>
#include <spg/sim/simObject.h>
#include <cmath>

namespace spg
{
namespace
{
auto l_stvkConstraint = [](const MembraneStvkEnergy *energy, const int i, const SimObject &obj, auto &dC) {
    const auto &x0p{obj.positions()[energy->stencils()[i][0]]};
    const auto &x1p{obj.positions()[energy->stencils()[i][1]]};
    const auto &x2p{obj.positions()[energy->stencils()[i][2]]};

    using RealT = std::decay_t<decltype(dC[0])>;
    const Vector3T<RealT> x0(RealT(x0p.x(), 0), RealT(x0p.y(), 1), RealT(x0p.z(), 2));
    const Vector3T<RealT> x1(RealT(x1p.x(), 3), RealT(x1p.y(), 4), RealT(x1p.z(), 5));
    const Vector3T<RealT> x2(RealT(x2p.x(), 6), RealT(x2p.y(), 7), RealT(x2p.z(), 8));
    // Compute deformation gradient
    const Vector3T<RealT> u = x1 - x0;
    const Vector3T<RealT> v = x2 - x0;
    const MatrixT<RealT, 3, 2> deformedMatrix{{u[0], v[0]}, {u[1], v[1]}, {u[2], v[2]}};
    const MatrixT<RealT, 3, 2> F = deformedMatrix * energy->inverseReferenceMats()[i];
    // Compute green strain
    // Note: Could be made more efficient by computing only required terms
    const MatrixT<RealT, 2, 2> E = 0.5 * (F.transpose() * F - Matrix2::Identity());

    dC = {E(0, 0), E(1, 1), 2 * E(0, 1)};
};
}

void MembraneStvkEnergy::addStencil(const std::array<int, s_stencilSize> &stencil, const Real young, const Real poisson)
{
    m_stencils.push_back(stencil);

    const Real youngCommon = young / (1 - poisson * poisson);
    const Matrix3 C{{youngCommon, youngCommon * poisson, 0},
                    {youngCommon * poisson, youngCommon, 0},
                    {0, 0, youngCommon * (1 - poisson) / 2}};
    m_modelStiffness.emplace_back(C);
    m_effectiveStiffness.emplace_back(C);
    m_modelCompliance.emplace_back(C.inverse());
    m_effectiveCompliance.emplace_back(C.inverse());
}

void MembraneStvkEnergy::preparePrecomputations(const SimObject &obj)
{
    const int nstencils{static_cast<int>(m_stencils.size())};
    m_inverseReferenceMat.resize(nstencils);
    for (int i = 0; i < nstencils; ++i) {
        const auto &x0{obj.positions0()[m_stencils[i][0]]};
        const auto &x1{obj.positions0()[m_stencils[i][1]]};
        const auto &x2{obj.positions0()[m_stencils[i][2]]};
        if (x0.z() != 0.0 || x1.z() != 0.0 || x2.z() != 0.0) {
            throw std::runtime_error("Rest positions of triangles for membrane StVK energy must have coordinate z=0");
        }
        // compute rest material matrix
        const Vector2 u = (x1 - x0).segment<2>(0);
        const Vector2 v = (x2 - x0).segment<2>(0);

        m_inverseReferenceMat[i].col(0) = u;
        m_inverseReferenceMat[i].col(1) = v;
        m_inverseReferenceMat[i] = m_inverseReferenceMat[i].inverse().eval();
        const auto area = (u.x() * v.y() - v.x() * u.y()) * 0.5;
        if (area <= 0) {
            throw std::runtime_error("Rest area of membrane StVK triangle is not positive");
        }
        m_effectiveCompliance[i] = m_modelCompliance[i] / area;
        m_effectiveStiffness[i] = m_modelStiffness[i] * area;
    }
    StencilBlockEnergy<3, 3>::preparePrecomputations(obj);
}

void MembraneStvkEnergy::dConstraints(int i, const SimObject &obj, ConstraintsAD1 &dC) const
{
    l_stvkConstraint(this, i, obj, dC);
}

void MembraneStvkEnergy::dConstraints(int i, const SimObject &obj, ConstraintsAD2 &dC) const
{
    l_stvkConstraint(this, i, obj, dC);
}
}  // namespace spg