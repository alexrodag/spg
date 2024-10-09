#include <spg/sim/energy/staticLoadEnergy.h>
#include <spg/sim/simObject.h>

namespace spg
{
namespace
{
auto l_loadConstraint =
    [](const StaticLoadEnergy * /*energy*/, const int /*i*/, const SimObject & /*obj*/, auto & /*dC*/) {
        // Note: I don't know if this energy can be implemented in constraint form. in PBD/XPBD these are typically
        // treated as external forces directly
        throw std::runtime_error("Constraint-based load energy not implemented");
    };

auto l_loadEnergy = [](const StaticLoadEnergy *energy, const int i, const SimObject &obj, auto &dE) {
    const auto &x0{obj.positions()[energy->stencils()[i][0]]};
    using ADouble = std::decay_t<decltype(dE)>;
    const Vector3T<ADouble> x(ADouble(x0.x(), 0), ADouble(x0.y(), 1), ADouble(x0.z(), 2));
    dE = energy->modelStiffness()[i][0] * -energy->loadDirections()[i].dot(x);
};
}  // namespace

void StaticLoadEnergy::dEnergy(const int i, const SimObject &obj, DScalarFirstD &dC) const
{
    l_loadEnergy(this, i, obj, dC);
}

void StaticLoadEnergy::dEnergy(const int i, const SimObject &obj, DScalarSecondD &dC) const
{
    l_loadEnergy(this, i, obj, dC);
}

void StaticLoadEnergy::dConstraints(const int i, const SimObject &obj, DConstraintsFirstD &dC) const
{
    l_loadConstraint(this, i, obj, dC);
}

void StaticLoadEnergy::dConstraints(const int i, const SimObject &obj, DConstraintsSecondD &dC) const
{
    l_loadConstraint(this, i, obj, dC);
}
}  // namespace spg