#include <spg/sim/solver/quasiStaticNewton.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/utils/timer.h>

#include <iostream>

namespace spg::solver
{
void QuasiStaticNewton::step()
{
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "QuasiStaticNewton step\n";
    }
    Timer timer;
    timer.start();

    const Real dt = m_dtStep / m_nsubsteps;
    const Real dtdt = dt * dt;
    // Compute total DOFs
    int accumulatedNDOF = 0;
    for (const auto &object : std::get<std::vector<SimObject>>(m_objects)) {
        accumulatedNDOF += object.nDOF();
    }
    const int totalNDOF{accumulatedNDOF};

    for (int s = 0; s < m_nsubsteps; ++s) {
        // Store state backup
        VectorX x0(totalNDOF);
        getSystemPositions(x0);

        // Set initial guess
        const VectorX xi = x0;
        // Note: No need to update objects positions since the guess is the current positions

        // Compute forces, mass matrix and stiffness matrix
        VectorX f(totalNDOF);
        SparseMatrix M(totalNDOF, totalNDOF);
        SparseMatrix K(totalNDOF, totalNDOF);
        getSystemForce(f);
        getSystemMassMatrix(M);
        getSystemStiffnessMatrix(K);

        // Create Linear problem left and right hand sides
        const SparseMatrix LHS = M - dtdt * K;
        const VectorX RHS = dtdt * f - M * (xi - x0);

        // Solve problem to obtain dx
        VectorX dx;
        solveLinearSystem(LHS, RHS, dx);

        // Update objects state
        const VectorX x = xi + dx;
        setObjectsPositions(x);
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Total step time: " << timer.getMilliseconds() << "ms\n\n";
    }
}
}  // namespace spg::solver