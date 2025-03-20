#pragma once

#include <spg/sim/energy/energy.h>

namespace spg
{
// FEM energy corresponding to Smith 2018, without the Log barrier term (i.e., corresponding to equation 11 of the
// paper), also in XPBD format (refs: "Stable Neo-Hookean Flesh Simulation", "Constraint-based Formulation of Stable
// Neo-Hookean Materials", "Parallel Block Neo-Hookean XPBD using Graph Clustering")
class StableNeoHookeanEnergy : public StencilBlockEnergy<4, 2>
{
public:
    StableNeoHookeanEnergy() { m_name = "StableNeoHookeanEnergy"; }
    void addStencil(const std::array<int, s_stencilSize> &stencil, Real young, Real poisson);
    virtual void preparePrecomputations(const ParticleGroup &pGroup);
    const std::vector<Matrix3> &inverseReferenceMats() const { return m_inverseReferenceMat; }
    const std::vector<Real> &hydrostaticCorrectionTerms() const { return m_hydrostaticCorrectionTerm; }

protected:
    virtual void dConstraints(int i, const ParticleGroup &pGroup, ConstraintsAD1 &dC) const;
    virtual void dConstraints(int i, const ParticleGroup &pGroup, ConstraintsAD2 &dC) const;
    virtual void dEnergy(int i, const ParticleGroup &pGroup, RealAD1 &dC) const;
    virtual void dEnergy(int i, const ParticleGroup &pGroup, RealAD2 &dC) const;
    std::vector<Matrix3> m_inverseReferenceMat;
    std::vector<Real> m_hydrostaticCorrectionTerm;
};

}  // namespace spg