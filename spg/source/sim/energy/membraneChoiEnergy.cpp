#include <spg/sim/energy/membraneChoiEnergy.h>
#include <spg/sim/simObject.h>
#include <cmath>

namespace spg
{
namespace
{
auto l_choiConstraint = [](const MembraneChoiEnergy *energy, const int i, const SimObject &obj, auto &dC) {
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
    const Real sqrtTwoInv = 1.0 / std::sqrt(2.0);
    dC = {F.col(0).norm() - 1,
          F.col(1).norm() - 1,
          ((F.col(0) + F.col(1)) * sqrtTwoInv).norm() - 1,
          ((F.col(0) - F.col(1)) * sqrtTwoInv).norm() - 1};
};
}  // namespace

void MembraneChoiEnergy::addStencil(const std::array<int, s_stencilSize> &stencil,
                                    const Real kWeft,
                                    const Real kWarp,
                                    const Real kBias)
{
    m_stencils.push_back(stencil);

    const Matrix4 C{{kWeft, 0, 0, 0}, {0, kWarp, 0, 0}, {0, 0, kBias, 0}, {0, 0, 0, kBias}};
    m_modelStiffness.emplace_back(C);
    m_effectiveStiffness.emplace_back(C);
    m_modelCompliance.emplace_back(C.inverse());
    m_effectiveCompliance.emplace_back(C.inverse());
}

void MembraneChoiEnergy::preparePrecomputations(const SimObject &obj)
{
    const int nstencils{static_cast<int>(m_stencils.size())};
    m_inverseReferenceMat.resize(nstencils);
    for (int i = 0; i < nstencils; ++i) {
        const auto &x0{obj.positions0()[m_stencils[i][0]]};
        const auto &x1{obj.positions0()[m_stencils[i][1]]};
        const auto &x2{obj.positions0()[m_stencils[i][2]]};
        if (x0.z() != 0.0 || x1.z() != 0.0 || x2.z() != 0.0) {
            throw std::runtime_error("Rest positions of triangles for membrane Choi energy must have coordinate z=0");
        }
        // compute rest material matrix
        const Vector2 u = (x1 - x0).segment<2>(0);
        const Vector2 v = (x2 - x0).segment<2>(0);

        m_inverseReferenceMat[i].col(0) = u;
        m_inverseReferenceMat[i].col(1) = v;
        m_inverseReferenceMat[i] = m_inverseReferenceMat[i].inverse().eval();
        const auto area = (u.x() * v.y() - v.x() * u.y()) * 0.5;
        if (area <= 0) {
            throw std::runtime_error("Rest area of membrane Choi triangle is not positive");
        }
        m_effectiveCompliance[i] = m_modelCompliance[i] / area;
        m_effectiveStiffness[i] = m_modelStiffness[i] * area;
    }
    StencilBlockEnergy<3, 4>::preparePrecomputations(obj);
}

void MembraneChoiEnergy::dConstraints(int i, const SimObject &obj, ConstraintsAD1 &dC) const
{
    l_choiConstraint(this, i, obj, dC);
}

void MembraneChoiEnergy::dConstraints(int i, const SimObject &obj, ConstraintsAD2 &dC) const
{
    l_choiConstraint(this, i, obj, dC);
}
}  // namespace spg