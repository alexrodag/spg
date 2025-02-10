#include <spg/sim/solver/baseSolver.h>
#include <spg/sim/simObject.h>

#include <iostream>

namespace spg::solver
{
void BaseSolver::addSimObject(const SimObject &object)
{
    m_simObjects.push_back(object);
}
void BaseSolver::addRigidBodyGroup(const RigidBodyGroup &rbGroup)
{
    m_rigidBodyGroups.push_back(rbGroup);
}

void BaseSolver::reset()
{
    for (auto &obj : m_simObjects) {
        obj.reset();
    }
    for (auto &rbg : m_rigidBodyGroups) {
        rbg.reset();
    }
}
}  // namespace spg::solver