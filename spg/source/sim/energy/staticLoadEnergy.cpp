#include <spg/sim/energy/staticLoadEnergy.h>
#include <spg/sim/simObject/particleGroup.h>

namespace spg
{
namespace
{
auto l_loadConstraint =
    [](const StaticLoadEnergy * /*energy*/, const int /*i*/, const ParticleGroup & /*pGroup*/, auto & /*dC*/) {
        // Note: I don't know if this energy can be implemented in constraint form. in PBD/XPBD these are typically
        // treated as external forces directly
        throw std::runtime_error("Constraint-based load energy not implemented");
    };

auto l_loadEnergy = [](const StaticLoadEnergy *energy, const int i, const ParticleGroup &pGroup, auto &dE) {
    const auto &x0{pGroup.positions()[energy->stencils()[i][0]]};
    using RealT = std::decay_t<decltype(dE)>;
    const Vector3T<RealT> x(RealT(x0.x(), 0), RealT(x0.y(), 1), RealT(x0.z(), 2));
    dE = energy->modelStiffness()[i][0] * -energy->loadDirections()[i].dot(x);
};
}  // namespace

void StaticLoadEnergy::dEnergy(const int i, const ParticleGroup &pGroup, RealAD1 &dC) const
{
    l_loadEnergy(this, i, pGroup, dC);
}

void StaticLoadEnergy::dEnergy(const int i, const ParticleGroup &pGroup, RealAD2 &dC) const
{
    l_loadEnergy(this, i, pGroup, dC);
}

void StaticLoadEnergy::dConstraints(const int i, const ParticleGroup &pGroup, ConstraintsAD1 &dC) const
{
    l_loadConstraint(this, i, pGroup, dC);
}

void StaticLoadEnergy::dConstraints(const int i, const ParticleGroup &pGroup, ConstraintsAD2 &dC) const
{
    l_loadConstraint(this, i, pGroup, dC);
}
}  // namespace spg