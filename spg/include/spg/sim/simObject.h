#pragma once

#include <spg/types.h>
#include <spg/geom/triangleMesh.h>

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
    const std::vector<Real> &invMasses() const { return m_w; }
    const std::vector<Real> &masses() const { return m_m; }
    std::vector<std::shared_ptr<EnergyT>> &energies() { return m_energies; }
    const std::vector<std::shared_ptr<EnergyT>> &energies() const { return m_energies; }
    int nDOF() const { return static_cast<int>(m_x.size()) * 3; }
    int nParticles() const { return static_cast<int>(m_x.size()); }

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

class RigidBodyGroup
{
public:
    static constexpr int s_nDOFs = 6;
    using EnergyT = Energy<RigidBodyGroup>;

    const std::vector<Vector3> &positions() const { return m_x; }
    std::vector<Vector3> &positions() { return m_x; }
    const std::vector<Vector3> &positions0() const { return m_x0; }
    const std::vector<Vector3> &orientations() const { return m_theta; }
    std::vector<Vector3> &orientations() { return m_theta; }
    const std::vector<Vector3> &orientations0() const { return m_theta0; }
    std::vector<Vector3> &velocities() { return m_v; }
    const std::vector<Vector3> &velocities() const { return m_v; }
    std::vector<Vector3> &omegas() { return m_omegas; }
    const std::vector<Vector3> &omegas() const { return m_omegas; }
    const std::vector<Real> &invMasses() const { return m_w; }
    const std::vector<Real> &masses() const { return m_m; }
    const std::vector<Matrix3> &invInertias() const { return m_inertiaInv; }
    const std::vector<Matrix3> &inertias() const { return m_inertia; }
    std::vector<Matrix3> &rotationMatrices() { return m_rotationMatrix; }
    const std::vector<Matrix3> &rotationMatrices() const { return m_rotationMatrix; }
    const std::vector<TriangleMesh> &visualMeshes() const { return m_visualMesh; }
    std::vector<std::shared_ptr<EnergyT>> &energies() { return m_energies; }
    const std::vector<std::shared_ptr<EnergyT>> &energies() const { return m_energies; }

    // TODO: May be needed to have two versions, one for generealized coordinates and other for velocities
    int nDOF() const { return static_cast<int>(m_x.size()) * 3 + m_theta.size() * 3; }
    int nBodies() const { return static_cast<int>(m_x.size()); }

    void addBody(const Vector3 &p,
                 const Vector3 &p0,
                 Real mass,
                 const Vector3 &theta,
                 const Vector3 &theta0,
                 const Matrix3 &localInertia,
                 const TriangleMesh &mesh);
    void addEnergy(std::shared_ptr<EnergyT> energy);
    void removeEnergy(std::shared_ptr<EnergyT> energy);
    void scaleMasses(Real scaleFactor);
    void setName(std::string name) { m_name = name; }
    void reset();

    void updateInertia(int i);
    void updateInertias();
    void updateRotationMatrix(int i);
    void updateRotationMatrices();
    void updateTheta(int i);
    void updateThetas();

protected:
    std::string m_name;  // Mostly for debug

    // TODO Review storage approach, see if VectorX or MatrixX3 make more sense
    // Center of mass position and velocity
    std::vector<Vector3> m_x;
    std::vector<Vector3> m_xInitial;
    std::vector<Vector3> m_x0;
    std::vector<Vector3> m_v;

    // Orientation and angular velocity
    std::vector<Vector3> m_theta;  // Orientation in axis angle.
    std::vector<Vector3> m_thetaInitial;
    std::vector<Vector3> m_theta0;
    std::vector<Vector3> m_omegas;          // Angular velocity
    std::vector<Matrix3> m_rotationMatrix;  // Rotation matrix corresponding to the orientation

    // Mass
    std::vector<Real> m_m;
    std::vector<Real> m_w;

    // Inertia
    std::vector<Matrix3> m_localInertia;
    std::vector<Matrix3> m_localInertiaInv;
    std::vector<Matrix3> m_inertia;
    std::vector<Matrix3> m_inertiaInv;

    // Visual mesh in rigid body space
    std::vector<TriangleMesh> m_visualMesh;

    // Energy data
    std::vector<std::shared_ptr<EnergyT>> m_energies;
};
}  // namespace spg