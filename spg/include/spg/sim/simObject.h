#pragma once

#include <spg/types.h>

#include <vector>
#include <memory>

namespace spg
{
class Energy;

// TODO templatize DOFs per particle
class SimObject
{
public:
    const std::vector<Vector3> &positions() const { return m_x; }
    std::vector<Vector3> &positions() { return m_x; }
    const std::vector<Vector3> &positions0() const { return m_x0; }
    std::vector<Vector3> &velocities() { return m_v; }
    const std::vector<Vector3> &velocities() const { return m_v; }
    const std::vector<Real> &invMasses() const { return m_w; }
    const std::vector<Real> &masses() const { return m_m; }
    std::vector<std::shared_ptr<Energy>> &energies() { return m_energies; }
    const std::vector<std::shared_ptr<Energy>> &energies() const { return m_energies; }
    int nDOF() const { return static_cast<int>(m_x.size()) * 3; }
    int nParticles() const { return static_cast<int>(m_x.size()); }

    void addParticle(const Vector3 &position, const Vector3 &position0, Real mass);
    void addEnergy(std::shared_ptr<Energy> energy);
    void removeEnergy(std::shared_ptr<Energy> energy);
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
    std::vector<std::shared_ptr<Energy>> m_energies;
};
}  // namespace spg