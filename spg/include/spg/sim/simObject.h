#pragma once

#include <spg/types.h>

#include <vector>
#include <memory>

namespace spg
{
template <class TSimObject>
class Energy;

// TODO templatize DOFs per particle
class SimObject  // TODO: Rename to Particles or something like that
{
public:
    static constexpr int s_nDOFs = 3;
    using EnergyT = Energy<SimObject>;

    const std::vector<Vector3> &positions() const { return m_x; }
    std::vector<Vector3> &positions() { return m_x; }
    const std::vector<Vector3> &positions0() const { return m_x0; }
    std::vector<Vector3> &velocities() { return m_v; }
    const std::vector<Vector3> &velocities() const { return m_v; }
    void getPositions(VectorX &pos, int offsetIndex) const;
    void getVelocities(VectorX &vel, int offsetIndex) const;
    void setPositions(const VectorX &pos, int offsetIndex);
    void setVelocities(const VectorX &vel, int offsetIndex);
    void updatePositionsFromDx(const VectorX &dx, int offsetIndex);
    void integrateVelocities(Real dt);
    void computeIntegratedVelocities(const VectorX &oldPos, int offsetIndex, Real invdt);
    const std::vector<Real> &invMasses() const { return m_w; }
    const std::vector<Real> &masses() const { return m_m; }
    std::vector<std::shared_ptr<EnergyT>> &energies() { return m_energies; }
    const std::vector<std::shared_ptr<EnergyT>> &energies() const { return m_energies; }
    int nDOF() const { return static_cast<int>(m_x.size()) * 3; }
    int nElements() const { return static_cast<int>(m_x.size()); }

    void addParticle(const Vector3 &position, const Vector3 &position0, Real mass);
    void addEnergy(std::shared_ptr<EnergyT> energy);
    void removeEnergy(std::shared_ptr<EnergyT> energy);
    void scaleMasses(Real scaleFactor);
    void setName(std::string name) { m_name = name; }
    void reset();

protected:
    std::string m_name;  // Mostly for debug
    // Particle data
    // TODO Review storage approach, see if VectorX or MatrixX3 make more sense
    std::vector<Vector3> m_x;
    std::vector<Vector3> m_xInitial;
    std::vector<Vector3> m_x0;
    std::vector<Vector3> m_v;
    std::vector<Real> m_m;
    std::vector<Real> m_w;

    // Energy data
    std::vector<std::shared_ptr<EnergyT>> m_energies;
};
}  // namespace spg