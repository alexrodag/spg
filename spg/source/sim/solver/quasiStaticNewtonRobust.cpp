#include <spg/sim/solver/quasiStaticNewtonRobust.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/simObject.h>
#include <spg/utils/timer.h>

#include <iostream>

namespace spg::solver
{
void QuasiStaticNewtonRobust::step()
{
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "QuasiStaticNewtonRobust step\n";
    }
    Timer timer;
    timer.start();

    Real stepResidual{0.};

    // Compute total DOFs
    int accumulatedNDOF = 0;
    for (const auto &object : m_simObjects) {
        accumulatedNDOF += object.nDOF();
    }
    const int totalNDOF{accumulatedNDOF};

    for (int s = 0; s < m_nsubsteps; ++s) {
        bool successfulStep{true};

        // Update current dt
        // TODO Increase when possible
        const Real maxDt = m_dtStep / m_nsubsteps;
        if (m_currentDt == -1. || m_currentDt > maxDt) {
            m_currentDt = maxDt;
        }
        const Real dt = m_currentDt;
        const Real dtdt = dt * dt;

        // Store state backup
        // TODO Check if worth storing them as members
        VectorX x0(totalNDOF);
        getSystemPositions(x0);

        int currentNewtonIterations = 0;
        // Set initial guess
        VectorX xi = x0;
        do {
            // Compute forces, mass matrix and stiffness matrix
            VectorX f(totalNDOF);
            SparseMatrix M(totalNDOF, totalNDOF);
            SparseMatrix K(totalNDOF, totalNDOF);
            getSystemForce(f);
            getSystemMassMatrix(M);
            getSystemStiffnessMatrix(K);

            // Create linear problem left and right hand sides
            const SparseMatrix LHS = M - dtdt * K;
            const VectorX RHS = dtdt * f - M * (xi - x0);
            const Real initialResidual = RHS.norm();

            // Solve problem to obtain dv
            VectorX dx;
            solveLinearSystem(LHS, RHS, dx);

            // Update objects using line search
            stepResidual = initialResidual;
            int lineSearchIteration = 0;
            Real alpha = 1.;
            while (stepResidual >= initialResidual && lineSearchIteration < m_maxLineSearchIterations) {
                // Update objects state
                const VectorX x = xi + dx * alpha;
                setObjectsPositions(x);
                // compute new RHS
                getSystemForce(f);
                stepResidual = (dtdt * f - M * (x - x0)).norm();
                /*std::cout << "  Line search with step " << alpha << ". Step residual: " << stepResidual
                          << ". Initial residual: " << initialResidual << "\n";*/
                if (stepResidual < initialResidual) {
                    xi = x;
                }
                alpha *= 0.5;
                lineSearchIteration++;
            }
            // Check if line search failed, which means that the newton iteration failed and the step failed
            if (stepResidual >= initialResidual && stepResidual > m_newtonThreshold) {
                successfulStep = false;
            }
            if (m_verbosity == Verbosity::Performance) {
                std::cout << " Newton iter " << currentNewtonIterations << " with dt=" << dt
                          << ", residual=" << stepResidual << " (threshold: " << m_newtonThreshold << ")\n";
            }
            currentNewtonIterations++;
        } while (stepResidual > m_newtonThreshold && successfulStep && currentNewtonIterations < m_maxNewtonIterations);
        // Check if newton step was successful
        if (stepResidual > m_newtonThreshold) {
            successfulStep = false;
        }
        if (!successfulStep) {
            setObjectsPositions(x0);
            m_currentDt *= m_timeSplittingFactor;
            m_consecutiveSuccessfulSteps = 0;
            if (m_verbosity == Verbosity::Performance) {
                std::cout << "  Step failed, new dt: " << m_currentDt << "\n";
            }
        } else {
            m_consecutiveSuccessfulSteps++;
            if (m_consecutiveSuccessfulSteps > m_requiredSuccessfulStepsToIncreaseDt && m_currentDt < maxDt) {
                m_currentDt /= m_timeSplittingFactor;
                if (m_verbosity == Verbosity::Performance) {
                    std::cout << "  Stable integration, new dt: " << m_currentDt << "\n";
                }
            }
        }
    }
    timer.stop();
    if (m_verbosity == Verbosity::Performance) {
        std::cout << "  Total step time: " << timer.getMilliseconds() << "ms\n\n";
    }
}

void QuasiStaticNewtonRobust::reset()
{
    ImplicitEulerBase::reset();
    m_currentDt = -1.;
    m_consecutiveSuccessfulSteps = 0;
}
}  // namespace spg::solver