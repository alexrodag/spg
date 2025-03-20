#include <spg/sim/solver/implicitEulerNewtonDv.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject/particleGroup.h>
#include <spg/sim/simObject/rigidBodyGroup.h>
#include <spg/utils/timer.h>
#include <spg/utils/functionalUtilities.h>

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

        // Compute mass matrix in initial state to prevent simulations with rigid bodies to explode
        SparseMatrix M(totalNDOF, totalNDOF);
        getSystemMassMatrix(M);

        // Set initial guess as inertial position
        integrateObjectsVelocities(dt);

        // Compute forces, and stiffness matrix
        VectorX f(totalNDOF);
        SparseMatrix K(totalNDOF, totalNDOF);
        getSystemForce(f);
        getSystemStiffnessMatrix(K);

        // Create Linear problem left and right hand sides
        const SparseMatrix LHS = M - (dt * dt) * K;
        const VectorX RHS = dt * f;

        // Solve problem to obtain dv
        VectorX dv;
        solveLinearSystem(LHS, RHS, dv);

        // Update objects state
        setObjectsPositions(x0);
        setObjectsVelocities(v0 + dv);
        integrateObjectsVelocities(dt);
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Total step time: " << timer.getMilliseconds() << "ms\n\n";
    }
}
}  // namespace spg::solver