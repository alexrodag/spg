#include <spg/sim/solver/bdf2.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/utils/timer.h>

#include <iostream>

namespace spg::solver
{
void BDF2::step()
{
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "BDF2 step\n";
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
        getSystemPositions(x0);
        getSystemVelocities(v0);

        // Init minus-one state if not initialized
        if (m_xMinusOne.size() != x0.size() || m_vMinusOne.size() != v0.size()) {
            m_xMinusOne = x0;
            m_vMinusOne = v0;
        }

        // Compute forces, mass matrix and stiffness matrix
        VectorX f(totalNDOF);
        SparseMatrix M(totalNDOF, totalNDOF);
        SparseMatrix K(totalNDOF, totalNDOF);
        getSystemForce(f);
        getSystemMassMatrix(M);
        getSystemStiffnessMatrix(K);

        // Create Linear problem left and right hand sides
        const SparseMatrix LHS = M - 4.0 / 9.0 * (dt * dt) * K;
        const VectorX RHS = 1.0 / 3.0 * M * (v0 - m_vMinusOne) +
                            2.0 / 3.0 * dt * (f + 1.0 / 3.0 * (K * (x0 - m_xMinusOne + 2 * dt * v0)));

        // Solve problem to obtain dv
        VectorX dv;
        solveLinearSystem(LHS, RHS, dv);

        // Update objects state
        const VectorX v = v0 + dv;
        const VectorX x = 4.0 / 3.0 * x0 - 1.0 / 3.0 * m_xMinusOne + 2.0 / 3.0 * v * dt;

        // Store previous state
        m_xMinusOne = x0;
        m_vMinusOne = v0;

        setObjectsPositions(x);
        setObjectsVelocities(v);
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Total step time: " << timer.getMilliseconds() << "ms\n\n";
    }
}

void BDF2::reset()
{
    ImplicitEulerBase::reset();
    m_vMinusOne.resize(0, 0);
    m_xMinusOne.resize(0, 0);
}
}  // namespace spg::solver