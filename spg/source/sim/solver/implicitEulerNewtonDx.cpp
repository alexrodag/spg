#include <spg/sim/solver/implicitEulerNewtonDx.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/sim/rigidBodyGroup.h>
#include <spg/utils/timer.h>
#include <spg/utils/functionalUtilities.h>

#include <iostream>

namespace spg::solver
{
void ImplicitEulerNewtonDx::step()
{
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "NewtonDx step\n";
    }
    Timer timer;
    timer.start();

    const Real dt = m_dtStep / m_nsubsteps;
    // Compute total DOFs
    int accumulatedNDOF = 0;
    apply_each(
        [&accumulatedNDOF](const auto &objs) {
            for (const auto &obj : objs) {
                accumulatedNDOF += obj.nDOF();
            }
        },
        m_objects);
    const int totalNDOF{accumulatedNDOF};

    for (int s = 0; s < m_nsubsteps; ++s) {
        // Store state backup
        VectorX x0(totalNDOF);
        VectorX v0(totalNDOF);
        getSystemPositions(x0);
        getSystemVelocities(v0);

        // Set initial guess as inertial position
        integrateObjectsPositions(dt);
        VectorX xi(totalNDOF);
        getSystemPositions(xi);

        // Compute forces, mass matrix and stiffness matrix
        VectorX f(totalNDOF);
        SparseMatrix M(totalNDOF, totalNDOF);
        SparseMatrix K(totalNDOF, totalNDOF);
        getSystemForce(f);
        getSystemMassMatrix(M);
        getSystemStiffnessMatrix(K);

        // Create Linear problem left and right hand sides
        const Real invdt = 1. / dt;
        const SparseMatrix LHS = (invdt * invdt) * M - K;
        const VectorX RHS = f;

        // Solve problem to obtain dx
        VectorX dx;
        solveLinearSystem(LHS, RHS, dx);

        // Update objects state
        setObjectsPositions(x0);
        updateObjectPositionsAndIntegrateVelocities(xi + dx, invdt);
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Total step time: " << timer.getMilliseconds() << "ms\n\n";
    }
}
}  // namespace spg::solver