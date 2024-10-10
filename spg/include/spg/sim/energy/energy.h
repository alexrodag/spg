#pragma once

#include <spg/types.h>
#include <spg/sim/simObject.h>

#include <TinyAD/Scalar.hh>
#include <TinyAD/Utils/HessianProjection.hh>
#include <tuple>
#include <variant>

namespace spg
{
class SimObject;

// TODO: Check if passing functors to some methods allow more flexibility in solver creation without having to add more
// functions to Energy
// TODO: Template the class to account for nDOFs

// TODO: Decide if the const parameters in the virtual functions (int i, real dt, etc) should be const in the
// declaration, since the definition is going to be also in this .h file and C++ suggests having equal signature for
// overriden functions
class Energy
{
public:
    virtual ~Energy() = default;
    virtual void projectPosition(int i, SimObject &obj, Real dt) const = 0;
    virtual void updateExplicitVelocities(int i, SimObject &obj, Real dt, bool atomicUpdate) const = 0;
    virtual std::vector<Triplet> negativeHessianTriplets(int i, const SimObject &obj, int offsetDOF) const = 0;
    virtual void accumulateForces(int i,
                                  const SimObject &obj,
                                  int offsetDOF,
                                  VectorX &forcesVector,
                                  bool atomicUpdate) const = 0;
    virtual void accumulateNegativeHessianDiagonalAtomic(int i,
                                                         const SimObject &obj,
                                                         int offsetDOF,
                                                         VectorX &hessianDiagonalVector) const = 0;
    // TODO: Change to generic template accumulators when Energy is templatized
    virtual void accumulateVertexForce(int stencilIdx,
                                       int vertexIdxInStencil,
                                       const SimObject &obj,
                                       Vector3 &force) const = 0;
    virtual void accumulateVertexHessian(int stencilIdx,
                                         int vertexIdxInStencil,
                                         const SimObject &obj,
                                         Matrix3 &hessian) const = 0;
    virtual int nStencils() const = 0;
    virtual int stencilSize() const = 0;
    virtual std::vector<int> flatStencils() const = 0;
    virtual int nDOFs() const = 0;
    virtual Real energy(int i, const SimObject &obj) const = 0;
    // Note: These two are generic, but slower since they rely on runtime sized VectorX and MatrixX. Avoid their use if
    // possible.
    virtual VectorX energyGradientGeneric(int i, const SimObject &obj) const = 0;
    virtual MatrixX energyHessianGeneric(int i, const SimObject &obj) const = 0;

    virtual void preparePrecomputations(const SimObject & /*obj*/)
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
template <int TstencilSize, int TnConstraints = 1, int TnDOFs = 3>
class StencilBlockEnergy : public Energy
{
protected:
    // Static constexpr and type definitions to expose template parameters to subclasses
    // TODO: Find more expressive names for the type definitions
    static constexpr int s_stencilSize{TstencilSize};
    static constexpr int s_DOFs{TnDOFs};
    // Autodiff scalar with first derivative computation
    using DScalarFirstD = TinyAD::Scalar<TstencilSize * TnDOFs, Real, false>;
    // Autodiff scalar with first and second derivative computation
    using DScalarSecondD = TinyAD::Scalar<TstencilSize * TnDOFs, Real, true>;
    // Autodiff vector of nCosntraints size with first derivative computation
    using DConstraintsFirstD = VectorT<DScalarFirstD, TnConstraints>;
    // Autodiff vector of nCosntraints size with first and second derivative computation
    using DConstraintsSecondD = VectorT<DScalarSecondD, TnConstraints>;
    using ConstraintsType = VectorT<Real, TnConstraints>;
    using ConstraintsGradType = MatrixT<Real, TstencilSize * TnDOFs, TnConstraints>;
    using EnergyGradType = VectorT<Real, TstencilSize * TnDOFs>;
    using EnergyHessType = MatrixT<Real, TstencilSize * TnDOFs, TstencilSize * TnDOFs>;
    using StiffnessType = MatrixT<Real, TnConstraints, TnConstraints>;

public:
    virtual int nStencils() const { return static_cast<int>(m_stencils.size()); }
    virtual int stencilSize() const { return TstencilSize; }
    virtual int nDOFs() const { return TnDOFs; }
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
    const std::vector<StiffnessType> &effectiveStiffness() const { return m_effectiveStiffness; }
    const std::vector<StiffnessType> &modelStiffness() const { return m_modelStiffness; }
    const std::vector<StiffnessType> &effectiveCompliance() const { return m_effectiveCompliance; }
    const std::vector<StiffnessType> &modelCompliance() const { return m_modelCompliance; }

    virtual void projectPosition(int i, SimObject &obj, Real dt) const
    {
        // XPBD projection method. Constraints are solved in coupled fashion (ref: "Parallel Block Neo-Hookean XPBD
        // using Graph Clustering"). When TnConstraints == 1, it reduces to base XPBD formulation
        constexpr Real smallEpsilon{static_cast<Real>(1e-10)};
        const auto [C, grad] = constraintsAndGradient(i, obj);
        VectorT<Real, TstencilSize * TnDOFs> W;
        for (int s = 0; s < TstencilSize; ++s) {
            const auto w = obj.invMasses()[m_stencils[i][s]];
            for (int m = 0; m < TnDOFs; ++m) {
                W(s * TnDOFs + m) = w;
            }
        }
        const ConstraintsGradType WgradC = W.asDiagonal() * grad;
        const auto A =
            grad.transpose() * WgradC + m_effectiveCompliance[i] / (dt * dt) + smallEpsilon * StiffnessType::Identity();
        // Note: simplified numerator (ref: "Small Steps in Physics Simulation"). Doing a single iteration per step
        // allows to make this simplification but more importantly, maintains the accuracy of the solver, making it
        // converge to the real solution (see "XPBD: Position-Based Simulation of Compliant Constrained Dynamics")
        const ConstraintsType deltaLambda = -A.inverse() * C;
        const auto deltaX = (WgradC * deltaLambda).eval();
        for (int s = 0; s < TstencilSize; ++s) {
            obj.positions()[m_stencils[i][s]] += deltaX.template segment<TnDOFs>(s * TnDOFs);
        }
    }

    virtual void updateExplicitVelocities(int i, SimObject &obj, Real dt, bool atomicUpdate) const
    {
        const auto grad = energyGradient(i, obj);
        VectorT<Real, TstencilSize * TnDOFs> W;
        for (int s = 0; s < TstencilSize; ++s) {
            const auto w = obj.invMasses()[m_stencils[i][s]];
            for (int m = 0; m < TnDOFs; ++m) {
                W(s * TnDOFs + m) = w;
            }
        }
        const EnergyGradType deltaV = W.asDiagonal() * grad * -dt;
        if (!atomicUpdate) {
            for (int s = 0; s < TstencilSize; ++s) {
                obj.velocities()[m_stencils[i][s]] += deltaV.template segment<TnDOFs>(s * TnDOFs);
            }
        } else {
            for (int s = 0; s < TstencilSize; ++s) {
                for (int d = 0; d < TnDOFs; ++d) {
#pragma omp atomic
                    obj.velocities()[m_stencils[i][s]](d) += deltaV(s * TnDOFs + d);
                }
            }
        }
    }

    virtual void accumulateForces(int i,
                                  const SimObject &obj,
                                  int offsetDOF,
                                  VectorX &forcesVector,
                                  bool atomicUpdate) const
    {
        const EnergyGradType f = -energyGradient(i, obj);
        const auto &stencil = m_stencils[i];
        if (!atomicUpdate) {
            for (int particleId = 0; particleId < TstencilSize; ++particleId) {
                forcesVector.segment<TnDOFs>(stencil[particleId] * TnDOFs + offsetDOF) +=
                    f.template segment<TnDOFs>(particleId * TnDOFs);
            }
        } else {
            for (int particleId = 0; particleId < TstencilSize; ++particleId) {
                for (int d = 0; d < TnDOFs; ++d) {
#pragma omp atomic
                    forcesVector(stencil[particleId] * TnDOFs + d + offsetDOF) += f(particleId * TnDOFs + d);
                }
            }
        }
    };

    virtual void accumulateNegativeHessianDiagonalAtomic(int i,
                                                         const SimObject &obj,
                                                         int offsetDOF,
                                                         VectorX &hessianDiagonalVector) const
    {
        const EnergyHessType h = -energyHessian(i, obj);
        const auto &stencil = m_stencils[i];
        for (int particleId = 0; particleId < TstencilSize; ++particleId) {
            for (int d = 0; d < TnDOFs; ++d) {
                const auto index = particleId * TnDOFs + d;
#pragma omp atomic
                hessianDiagonalVector(stencil[particleId] * TnDOFs + d + offsetDOF) += h(index, index);
            }
        }
    };

    virtual std::vector<Triplet> negativeHessianTriplets(int i, const SimObject &obj, int offsetDOF) const
    {
        std::vector<Triplet> triplets;
        triplets.reserve(TstencilSize * TnDOFs * TstencilSize * TnDOFs);
        const EnergyHessType h = -energyHessian(i, obj);
        const auto &stencil = m_stencils[i];
        for (int rowStencil = 0; rowStencil < TstencilSize; ++rowStencil) {
            for (int colStencil = 0; colStencil < TstencilSize; ++colStencil) {
                for (int rowDOF = 0; rowDOF < TnDOFs; ++rowDOF) {
                    for (int colDOF = 0; colDOF < TnDOFs; ++colDOF) {
                        triplets.emplace_back(stencil[rowStencil] * TnDOFs + rowDOF + offsetDOF,
                                              stencil[colStencil] * TnDOFs + colDOF + offsetDOF,
                                              h(rowStencil * TnDOFs + rowDOF, colStencil * TnDOFs + colDOF));
                    }
                }
            }
        }
        return triplets;
    }

    // TODO: Check if possible to do autodiff only on required vertex in a generic way
    virtual void accumulateVertexForce(int stencilIdx,
                                       int vertexIdxInStencil,
                                       const SimObject &obj,
                                       Vector3 &force) const
    {
        static_assert(TnDOFs == 3, "accumulateVertexForce called with DoFs different than 3");
        const EnergyGradType f = -energyGradient(stencilIdx, obj);
        force += f.template segment<TnDOFs>(vertexIdxInStencil * TnDOFs);
    }

    // TODO: The current implementation computes the full hessian and then uses only that of one vertex. Check if
    // possible to do autodiff only on required vertex in a generic way
    virtual void accumulateVertexHessian(int stencilIdx,
                                         int vertexIdxInStencil,
                                         const SimObject &obj,
                                         Matrix3 &hessian) const
    {
        static_assert(TnDOFs == 3, "accumulateVertexForce called with DoFs different than 3");
        const EnergyHessType h = energyHessian(stencilIdx, obj);
        hessian += h.template block<TnDOFs, TnDOFs>(vertexIdxInStencil * TnDOFs, vertexIdxInStencil * TnDOFs);
    }

    // Note These core methods default to use the automatic differentiable versions but can be overloaded with more
    // efficient custom implementations. There are some examples in the implemented energies (e.g., SpringEnergy)
    virtual Real energy(int i, const SimObject &obj) const
    {
        DScalarFirstD dE;
        dEnergy(i, obj, dE);
        return dE.val;
    }
    virtual EnergyGradType energyGradient(int i, const SimObject &obj) const { return dEnergyGradient(i, obj); }

    virtual EnergyHessType energyHessian(int i, const SimObject &obj) const { return dEnergyHessian(i, obj); }

    virtual ConstraintsType constraints(int i, const SimObject &obj) const
    {
        DConstraintsFirstD dC;
        dConstraints(i, obj, dC);
        ConstraintsType C;
        // TODO: Make this conversion faster.
        for (int j = 0; j < TnConstraints; ++j) {
            C[j] = dC[j].val;
        }
        return C;
    }

    virtual ConstraintsGradType constraintsGradient(int i, const SimObject &obj) const
    {
        return std::get<1>(dConstraintsAndGradient(i, obj));
    }

    virtual std::tuple<ConstraintsType, ConstraintsGradType> constraintsAndGradient(int i, SimObject &obj) const
    {
        return dConstraintsAndGradient(i, obj);
    }

    virtual VectorX energyGradientGeneric(int i, const SimObject &obj) const { return energyGradient(i, obj); }

    virtual MatrixX energyHessianGeneric(int i, const SimObject &obj) const { return energyHessian(i, obj); }

protected:
    StencilBlockEnergy() = default;  // Protected constructor to prevent class instantiation

    // If dConstraints is provided, the energy and derivatives will be automatically computed. It is not mandatory to
    // implement if the energy version is provided and no constraint-based solver are to be used
    // Note: I use [[maybe_unused]] because of the combined declaration and definition in the header file that impedes
    // hiding the argument names to prevent the "unused variable" compile error
    virtual void dConstraints([[maybe_unused]] int i,
                              [[maybe_unused]] const SimObject &obj,
                              [[maybe_unused]] DConstraintsFirstD &dC) const
    {
        throw std::runtime_error("dConstraints not implemented for " + m_name + " energy.");
    }
    virtual void dConstraints([[maybe_unused]] int i,
                              [[maybe_unused]] const SimObject &obj,
                              [[maybe_unused]] DConstraintsSecondD &dC) const
    {
        throw std::runtime_error("dConstraints not implemented for " + m_name + " energy.");
    }
    virtual std::tuple<ConstraintsType, ConstraintsGradType> dConstraintsAndGradient(int i,
                                                                                     const SimObject &obj) const final
    {
        DConstraintsFirstD dC;
        dConstraints(i, obj, dC);
        ConstraintsType C;
        // TODO: Make this conversion faster.
        for (int j = 0; j < TnConstraints; ++j) {
            C[j] = dC[j].val;
        }
        ConstraintsGradType gradC;
        // TODO: Make this conversion faster.
        for (int j = 0; j < TnConstraints; ++j) {
            gradC.col(j) = dC[j].grad;
        }
        return {C, gradC};
    }

    // dEnergy defaults to use the constraint-based approach. It can be overloaded to define a direct differentiable
    // version to be more efficient, or in cases where a constraint-based formulation is not available
    virtual void dEnergy(int i, const SimObject &obj, DScalarFirstD &dE) const
    {
        DConstraintsFirstD dC;
        dConstraints(i, obj, dC);
        dE = 0.5 * dC.dot(this->m_effectiveStiffness[i] * dC);
    }
    virtual void dEnergy(int i, const SimObject &obj, DScalarSecondD &dE) const
    {
        DConstraintsSecondD dC;
        dConstraints(i, obj, dC);
        dE = 0.5 * dC.dot(this->m_effectiveStiffness[i] * dC);
    }
    virtual EnergyGradType dEnergyGradient(int i, const SimObject &obj) const final
    {
        DScalarFirstD dE;
        dEnergy(i, obj, dE);
        return dE.grad;
    }
    virtual EnergyHessType dEnergyHessian(int i, const SimObject &obj, bool projectPositiveDefinite = false) const final
    {
        DScalarSecondD dE;
        dEnergy(i, obj, dE);
        if (projectPositiveDefinite) {
            TinyAD::project_positive_definite(dE.Hess, Real(TinyAD::default_hessian_projection_eps));
        }
        return dE.Hess;
    }

    std::vector<std::array<int, TstencilSize>> m_stencils;
    std::vector<StiffnessType> m_effectiveStiffness;   // Practical stiffness including integration along element
    std::vector<StiffnessType> m_modelStiffness;       // Stiffness matrix of the energy, inverse of the base compliance
    std::vector<StiffnessType> m_effectiveCompliance;  // Practical compliance including integration along element
    std::vector<StiffnessType> m_modelCompliance;      // Base compliance, equivalent to the inverse of the stiffness
};
}  // namespace spg