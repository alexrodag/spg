#include <spg/sim/solver/implicitEulerBaraffWitkin.h>
#include <spg/sim/simObject/particleGroup.h>
#include <spg/sim/simObject/rigidBodyGroup.h>
#include <spg/utils/timer.h>
#include <spg/utils/functionalUtilities.h>

#include <iostream>

namespace spg::solver
{
void ImplicitEulerBaraffWitkin::step()
{
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "BaraffWitkin step\n";
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
        // Store state velocities
        VectorX v0(totalNDOF);
        getSystemVelocities(v0);

        // Compute forces, mass matrix and stiffness matrix
        VectorX f(totalNDOF);
        SparseMatrix M(totalNDOF, totalNDOF);
        SparseMatrix K(totalNDOF, totalNDOF);
        getSystemForce(f);
        getSystemMassMatrix(M);
        getSystemStiffnessMatrix(K);

        // Create Linear problem left and right hand sides
        const SparseMatrix LHS = M - (dt * dt) * K;
        const VectorX RHS = dt * (f + dt * (K * v0));

        // Solve problem to obtain dv
        VectorX dv;
        solveLinearSystem(LHS, RHS, dv);

        // Update objects state
        setObjectsVelocities(v0 + dv);
        integrateObjectsVelocities(dt);
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Total step time: " << timer.getMilliseconds() << "ms\n\n";
    }
}
}  // namespace spg::solver