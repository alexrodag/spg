#include <spg/sim/solver/baseSolver.h>
#include <spg/sim/simObject/particleGroup.h>
#include <spg/sim/simObject/rigidBodyGroup.h>
#include <spg/utils/functionalUtilities.h>

#include <iostream>

namespace spg::solver
{
template <typename TObj>
void BaseSolver::addObject(const TObj &object)
{
    std::get<std::vector<TObj>>(m_objects).push_back(object);
}

int BaseSolver::numSimObjects()
{
    int nObjects = 0;
    apply_each([&nObjects](auto &objs) { nObjects += objs.size(); }, m_objects);
    return nObjects;
}

void BaseSolver::reset()
{
    apply_each(
        [](auto &objs) {
            for (auto &obj : objs) {
                obj.reset();
            }
        },
        m_objects);
}

template void BaseSolver::addObject(const ParticleGroup &object);
template void BaseSolver::addObject(const RigidBodyGroup &object);
}  // namespace spg::solver