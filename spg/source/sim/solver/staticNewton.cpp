#include <spg/sim/solver/staticNewton.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/sim/rigidBodyGroup.h>
#include <spg/utils/timer.h>

#include <iostream>

namespace spg::solver
{
void StaticNewton::step()
{
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "StaticNewton step\n";
    }
    Timer timer;
    timer.start();

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
        setObjectsPositions(xi);

        // Compute forces, mass matrix and stiffness matrix
        VectorX f(totalNDOF);
        SparseMatrix K(totalNDOF, totalNDOF);
        getSystemForce(f);
        getSystemStiffnessMatrix(K);

        // Create Linear problem left and right hand sides
        const SparseMatrix LHS = K;
        const VectorX RHS = -f;

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