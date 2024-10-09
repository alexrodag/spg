#include <spg/sim/solver/baseSolver.h>
#include <spg/sim/simObject.h>

#include <iostream>

namespace spg::solver
{
void BaseSolver::addSimObject(const SimObject &object)
{
    m_simObjects.push_back(object);
}

void BaseSolver::reset()
{
    for (auto &obj : m_simObjects) {
        obj.reset();
    }
}
}  // namespace spg::solver