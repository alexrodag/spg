#pragma once

#include <spg/types.h>
#include <spg/sim/simObject.h>

#include <TinyAD/Scalar.hh>
#include <TinyAD/Utils/HessianProjection.hh>
#include <tuple>

namespace spg
{

// Base Energy class providing a common interface
// TODO: Check if passing functors to some methods allow more flexibility in solver creation without having to add more
// functions to Energy
// TODO: Template the class to account for nDOFs
// TODO: Decide if the const parameters in the virtual functions (int i, real dt, etc) should be const in the
// declaration, since the definition is going to be also in this .h file and C++ suggests having equal signature for
// overriden functions
template <class TSimObject>
class Energy
{
public:
    static constexpr int s_nDOFs = TSimObject::s_nDOFs;
    virtual ~Energy() = default;
    virtual void projectPosition(int i, TSimObject &obj, Real dt) const = 0;
    virtual void updateExplicitVelocities(int i, TSimObject &obj, Real dt, bool atomicUpdate) const = 0;
    virtual std::vector<Triplet> negativeHessianTriplets(int i, const TSimObject &obj, int offsetDOF) const = 0;
    virtual void accumulateForces(int i,
                                  const TSimObject &obj,
                                  int offsetDOF,
                                  VectorX &forcesVector,
                                  bool atomicUpdate) const = 0;
    // TODO: Change to generic template accumulators when Energy is templatized
    virtual void accumulateVertexForce(int stencilIdx,
                                       int vertexIdxInStencil,
                                       const TSimObject &obj,
                                       Vector<s_nDOFs> &force) const = 0;
    virtual void accumulateVertexHessian(int stencilIdx,
                                         int vertexIdxInStencil,
                                         const TSimObject &obj,
                                         Matrix<s_nDOFs> &hessian) const = 0;
    virtual int nStencils() const = 0;
    virtual int stencilSize() const = 0;
    virtual std::vector<int> flatStencils() const = 0;
    virtual int nDOFs() const = 0;
    virtual Real energy(int i, const TSimObject &obj) const = 0;
    // Note: These two are generic, but slower since they rely on runtime sized VectorX and MatrixX. Avoid their use if
    // possible.
    virtual VectorX energyGradientGeneric(int i, const TSimObject &obj) const = 0;
    virtual MatrixX energyHessianGeneric(int i, const TSimObject &obj) const = 0;

    virtual void preparePrecomputations(const TSimObject & /*obj*/)
    {
        // Note: Currently I assume a single preparation is accepted. Could change in the future.
        if (m_prepared) {
            throw std::runtime_error("Energy was already prepared");
        }
        m_prepared = true;
    }

protected:
    bool m_prepared{false};
    std::string m_name;  // Mostly for debug
};

// This base class implements the core logic to do automatic differentiation of constraint- and energy-based energies.
// By providing just the constraint (or energy) implementation, all the energy and corresponding derivatives are
// automatically provided, as well as solver-specific functions. Through overloading, it is also possible to provide
// more efficient implementations
template <int TstencilSize, int TnConstraints = 1, class TSimObject = SimObject>
class StencilBlockEnergy : public Energy<TSimObject>
{
protected:
    // Static constexpr and type definitions to expose template parameters to subclasses
    static constexpr int s_nDOFs = Energy::s_nDOFs;
    static constexpr int s_stencilSize{TstencilSize};

    // Autodiff scalar with first derivative computation
    using RealAD1 = TinyAD::Scalar<TstencilSize * s_nDOFs, Real, false>;
    // Autodiff scalar with first and second derivative computation
    using RealAD2 = TinyAD::Scalar<TstencilSize * s_nDOFs, Real, true>;
    // Autodiff vector of nCosntraints size with first derivative computation
    using ConstraintsAD1 = VectorT<RealAD1, TnConstraints>;
    // Autodiff vector of nCosntraints size with first and second derivative computation
    using ConstraintsAD2 = VectorT<RealAD2, TnConstraints>;

    using Constraints = VectorT<Real, TnConstraints>;
    using ConstraintsGrad = MatrixT<Real, TstencilSize * s_nDOFs, TnConstraints>;
    using EnergyGrad = VectorT<Real, TstencilSize * s_nDOFs>;
    using EnergyHess = MatrixT<Real, TstencilSize * s_nDOFs, TstencilSize * s_nDOFs>;
    using StiffnessMat = MatrixT<Real, TnConstraints, TnConstraints>;

public:
    virtual int nStencils() const { return static_cast<int>(m_stencils.size()); }
    virtual int stencilSize() const { return TstencilSize; }
    virtual int nDOFs() const { return s_nDOFs; }
    virtual std::vector<int> flatStencils() const
    {
        std::vector<int> flattened;
        flattened.reserve(m_stencils.size() * TstencilSize);
        for (const auto &stencil : m_stencils) {
            flattened.insert(flattened.end(), stencil.begin(), stencil.end());
        }
        return flattened;
    }
    const std::vector<std::array<int, TstencilSize>> &stencils() const { return m_stencils; }
    const std::vector<StiffnessMat> &effectiveStiffness() const { return m_effectiveStiffness; }
    const std::vector<StiffnessMat> &modelStiffness() const { return m_modelStiffness; }
    const std::vector<StiffnessMat> &effectiveCompliance() const { return m_effectiveCompliance; }
    const std::vector<StiffnessMat> &modelCompliance() const { return m_modelCompliance; }

    virtual void projectPosition(int i, TSimObject &obj, Real dt) const
    {
        if constexpr (TSimObject::s_nDOFs == 3) {
            // XPBD projection method. Constraints are solved in coupled fashion (ref: "Parallel Block Neo-Hookean XPBD
            // using Graph Clustering"). When TnConstraints == 1, it reduces to base XPBD formulation
            constexpr Real smallEpsilon{static_cast<Real>(1e-10)};
            const auto [C, grad] = constraintsAndGradient(i, obj);
            VectorT<Real, TstencilSize * s_nDOFs> W;
            for (int s = 0; s < TstencilSize; ++s) {
                const auto w = obj.invMasses()[m_stencils[i][s]];
                for (int m = 0; m < s_nDOFs; ++m) {
                    W(s * s_nDOFs + m) = w;
                }
            }
            const ConstraintsGrad WgradC = W.asDiagonal() * grad;
            const auto A = grad.transpose() * WgradC + m_effectiveCompliance[i] / (dt * dt) +
                           smallEpsilon * StiffnessMat::Identity();
            // Note: simplified numerator (ref: "Small Steps in Physics Simulation"). Doing a single iteration per step
            // allows to make this simplification but more importantly, maintains the accuracy of the solver, making it
            // converge to the real solution (see "XPBD: Position-Based Simulation of Compliant Constrained Dynamics")
            const Constraints deltaLambda = -A.inverse() * C;
            const auto deltaX = (WgradC * deltaLambda).eval();
            for (int s = 0; s < TstencilSize; ++s) {
                obj.positions()[m_stencils[i][s]] += deltaX.template segment<s_nDOFs>(s * s_nDOFs);
            }
        }
    }

    virtual void updateExplicitVelocities(int i, TSimObject &obj, Real dt, bool atomicUpdate) const
    {
        if constexpr (TSimObject::s_nDOFs == 3) {
            const auto grad = energyGradient(i, obj);
            VectorT<Real, TstencilSize * s_nDOFs> W;
            for (int s = 0; s < TstencilSize; ++s) {
                const auto w = obj.invMasses()[m_stencils[i][s]];
                for (int m = 0; m < s_nDOFs; ++m) {
                    W(s * s_nDOFs + m) = w;
                }
            }
            const EnergyGrad deltaV = W.asDiagonal() * grad * -dt;
            if (!atomicUpdate) {
                for (int s = 0; s < TstencilSize; ++s) {
                    obj.velocities()[m_stencils[i][s]] += deltaV.template segment<s_nDOFs>(s * s_nDOFs);
                }
            } else {
                for (int s = 0; s < TstencilSize; ++s) {
                    for (int d = 0; d < s_nDOFs; ++d) {
#pragma omp atomic
                        obj.velocities()[m_stencils[i][s]](d) += deltaV(s * s_nDOFs + d);
                    }
                }
            }
        }
    }

    virtual void accumulateForces(int i,
                                  const TSimObject &obj,
                                  int offsetDOF,
                                  VectorX &forcesVector,
                                  bool atomicUpdate) const
    {
        const EnergyGrad f = -energyGradient(i, obj);
        const auto &stencil = m_stencils[i];
        if (!atomicUpdate) {
            for (int particleId = 0; particleId < TstencilSize; ++particleId) {
                forcesVector.segment<s_nDOFs>(stencil[particleId] * s_nDOFs + offsetDOF) +=
                    f.template segment<s_nDOFs>(particleId * s_nDOFs);
            }
        } else {
            for (int particleId = 0; particleId < TstencilSize; ++particleId) {
                for (int d = 0; d < s_nDOFs; ++d) {
#pragma omp atomic
                    forcesVector(stencil[particleId] * s_nDOFs + d + offsetDOF) += f(particleId * s_nDOFs + d);
                }
            }
        }
    };

    virtual std::vector<Triplet> negativeHessianTriplets(int i, const TSimObject &obj, int offsetDOF) const
    {
        std::vector<Triplet> triplets;
        triplets.reserve(TstencilSize * s_nDOFs * TstencilSize * s_nDOFs);
        const EnergyHess h = -energyHessian(i, obj);
        const auto &stencil = m_stencils[i];
        for (int rowStencil = 0; rowStencil < TstencilSize; ++rowStencil) {
            for (int colStencil = 0; colStencil < TstencilSize; ++colStencil) {
                for (int rowDOF = 0; rowDOF < s_nDOFs; ++rowDOF) {
                    for (int colDOF = 0; colDOF < s_nDOFs; ++colDOF) {
                        triplets.emplace_back(stencil[rowStencil] * s_nDOFs + rowDOF + offsetDOF,
                                              stencil[colStencil] * s_nDOFs + colDOF + offsetDOF,
                                              h(rowStencil * s_nDOFs + rowDOF, colStencil * s_nDOFs + colDOF));
                    }
                }
            }
        }
        return triplets;
    }

    // TODO: Check if possible to do autodiff only on required vertex in a generic way
    virtual void accumulateVertexForce(int stencilIdx,
                                       int vertexIdxInStencil,
                                       const TSimObject &obj,
                                       Vector<s_nDOFs> &force) const
    {
        const EnergyGrad f = -energyGradient(stencilIdx, obj);
        force += f.template segment<s_nDOFs>(vertexIdxInStencil * s_nDOFs);
    }

    // TODO: The current implementation computes the full hessian and then uses only that of one vertex. Check if
    // possible to do autodiff only on required vertex in a generic way
    virtual void accumulateVertexHessian(int stencilIdx,
                                         int vertexIdxInStencil,
                                         const TSimObject &obj,
                                         Matrix<s_nDOFs> &hessian) const
    {
        const EnergyHess h = energyHessian(stencilIdx, obj);
        hessian += h.template block<s_nDOFs, s_nDOFs>(vertexIdxInStencil * s_nDOFs, vertexIdxInStencil * s_nDOFs);
    }

    // Note These core methods default to use the automatic differentiable versions but can be overloaded with more
    // efficient custom implementations. There are some examples in the implemented energies (e.g., SpringEnergy)
    virtual Real energy(int i, const TSimObject &obj) const
    {
        RealAD1 dE;
        dEnergy(i, obj, dE);
        return dE.val;
    }
    virtual EnergyGrad energyGradient(int i, const TSimObject &obj) const { return dEnergyGradient(i, obj); }

    virtual EnergyHess energyHessian(int i, const TSimObject &obj) const { return dEnergyHessian(i, obj); }

    virtual Constraints constraints(int i, const TSimObject &obj) const
    {
        ConstraintsAD1 dC;
        dConstraints(i, obj, dC);
        Constraints C;
        // TODO: Make this conversion faster.
        for (int j = 0; j < TnConstraints; ++j) {
            C[j] = dC[j].val;
        }
        return C;
    }

    virtual ConstraintsGrad constraintsGradient(int i, const TSimObject &obj) const
    {
        return std::get<1>(dConstraintsAndGradient(i, obj));
    }

    virtual std::tuple<Constraints, ConstraintsGrad> constraintsAndGradient(int i, TSimObject &obj) const
    {
        return dConstraintsAndGradient(i, obj);
    }

    virtual VectorX energyGradientGeneric(int i, const TSimObject &obj) const { return energyGradient(i, obj); }

    virtual MatrixX energyHessianGeneric(int i, const TSimObject &obj) const { return energyHessian(i, obj); }

protected:
    StencilBlockEnergy() = default;  // Protected constructor to prevent class instantiation

    // If dConstraints is provided, the energy and derivatives will be automatically computed. It is not mandatory to
    // implement if the energy version is provided and no constraint-based solver are to be used
    // Note: I use [[maybe_unused]] because of the combined declaration and definition in the header file that impedes
    // hiding the argument names to prevent the "unused variable" compile error
    virtual void dConstraints([[maybe_unused]] int i,
                              [[maybe_unused]] const TSimObject &obj,
                              [[maybe_unused]] ConstraintsAD1 &dC) const
    {
        throw std::runtime_error("dConstraints not implemented for " + m_name + " energy.");
    }
    virtual void dConstraints([[maybe_unused]] int i,
                              [[maybe_unused]] const TSimObject &obj,
                              [[maybe_unused]] ConstraintsAD2 &dC) const
    {
        throw std::runtime_error("dConstraints not implemented for " + m_name + " energy.");
    }
    virtual std::tuple<Constraints, ConstraintsGrad> dConstraintsAndGradient(int i, const TSimObject &obj) const final
    {
        ConstraintsAD1 dC;
        dConstraints(i, obj, dC);
        Constraints C;
        // TODO: Make this conversion faster.
        for (int j = 0; j < TnConstraints; ++j) {
            C[j] = dC[j].val;
        }
        ConstraintsGrad gradC;
        // TODO: Make this conversion faster.
        for (int j = 0; j < TnConstraints; ++j) {
            gradC.col(j) = dC[j].grad;
        }
        return {C, gradC};
    }

    // dEnergy defaults to use the constraint-based approach. It can be overloaded to define a direct differentiable
    // version to be more efficient, or in cases where a constraint-based formulation is not available
    virtual void dEnergy(int i, const TSimObject &obj, RealAD1 &dE) const
    {
        ConstraintsAD1 dC;
        dConstraints(i, obj, dC);
        dE = 0.5 * dC.dot(this->m_effectiveStiffness[i] * dC);
    }
    virtual void dEnergy(int i, const TSimObject &obj, RealAD2 &dE) const
    {
        ConstraintsAD2 dC;
        dConstraints(i, obj, dC);
        dE = 0.5 * dC.dot(this->m_effectiveStiffness[i] * dC);
    }
    virtual EnergyGrad dEnergyGradient(int i, const TSimObject &obj) const final
    {
        RealAD1 dE;
        dEnergy(i, obj, dE);
        return dE.grad;
    }
    virtual EnergyHess dEnergyHessian(int i, const TSimObject &obj, bool projectPositiveDefinite = false) const final
    {
        RealAD2 dE;
        dEnergy(i, obj, dE);
        if (projectPositiveDefinite) {
            TinyAD::project_positive_definite(dE.Hess, Real(TinyAD::default_hessian_projection_eps));
        }
        return dE.Hess;
    }

    std::vector<std::array<int, TstencilSize>> m_stencils;
    std::vector<StiffnessMat> m_effectiveStiffness;   // Practical stiffness including integration along element
    std::vector<StiffnessMat> m_modelStiffness;       // Stiffness matrix of the energy, inverse of the base compliance
    std::vector<StiffnessMat> m_effectiveCompliance;  // Practical compliance including integration along element
    std::vector<StiffnessMat> m_modelCompliance;      // Base compliance, equivalent to the inverse of the stiffness
};
}  // namespace spg