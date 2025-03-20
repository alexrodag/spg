#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
// Bending energy corresponding to Bergou 2006
// (ref: "A Quadratic Bending Model for Inextensible Surfaces")
class QuadraticBendingEnergy : public StencilBlockEnergy<4, 3>
{
public:
    QuadraticBendingEnergy() { m_name = "QuadraticBendingEnergy"; }
    void addStencil(const std::array<int, s_stencilSize> &stencil, Real restTheta, Real stiffness);
    virtual void preparePrecomputations(const ParticleGroup &obj);
    const std::vector<Matrix<3, 12>> &laplacians() const { return m_laplacians; }

protected:
    virtual void dConstraints(int i, const ParticleGroup &obj, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const ParticleGroup &obj, ConstraintsAD2 &dC) const;
    std::vector<Matrix<3, 12>> m_laplacians;
    std::vector<Real> m_stiffnessFactor;
};

}  // namespace spg