#include <spg/sim/solver/implicitEulerNewtonDv.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/utils/timer.h>

#include <iostream>

namespace spg::solver
{
void ImplicitEulerNewtonDv::step()
{
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "NewtonDv step\n";
    }
    Timer timer;
    timer.start();

    const Real dt = m_dtStep / m_nsubsteps;
    // Compute total DOFs
    int accumulatedNDOF = 0;
    for (const auto &object : m_simObjects) {
        accumulatedNDOF += object.nDOF();
    }
    const int totalNDOF{accumulatedNDOF};

    for (int s = 0; s < m_nsubsteps; ++s) {
        // Store state backup
        VectorX x0(totalNDOF);
        VectorX v0(totalNDOF);
        getSystemState(x0, v0);

        // Set initial guess
        const VectorX xi = x0 + dt * v0;
        const VectorX vi = v0;
        setObjectsPositions(xi);

        // Compute forces, mass matrix and stiffness matrix
        VectorX f(totalNDOF);
        SparseMatrix M(totalNDOF, totalNDOF);
        SparseMatrix K(totalNDOF, totalNDOF);
        getSystemForce(f);
        getSystemMassMatrix(M);
        getSystemStiffnessMatrix(K);

        // Create Linear problem left and right hand sides
        const SparseMatrix LHS = M - (dt * dt) * K;
        const VectorX RHS = dt * f;

        // Solve problem to obtain dv
        VectorX dv;
        solveLinearSystem(LHS, RHS, dv);

        // Update objects state
        const VectorX v = v0 + dv;
        const VectorX x = x0 + v * dt;
        setObjectsState(x, v);
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Total step time: " << timer.getMilliseconds() << "ms\n\n";
    }
}
}  // namespace spg::solver