#include <spg/sim/simObject.h>
#include <spg/sim/energy/energy.h>

namespace
{
spg::Matrix3 axisAngleToRotMatrix(const spg::Vector3 &theta)
{
    const auto angle = theta.norm();
    return angle != 0 ? Eigen::AngleAxis<spg::Real>(angle, theta / angle).toRotationMatrix() : spg::Matrix3::Identity();
}

spg::Vector3 rotMatrixToAxisAngle(const spg::Matrix3 &R)
{
    Eigen::AngleAxis<spg::Real> axisAngleOrientation(R);
    return axisAngleOrientation.angle() != 0 ? (axisAngleOrientation.axis() / axisAngleOrientation.angle()).eval()
                                             : spg::Vector3::Zero();
}
}  // namespace
namespace spg
{
void SimObject::addParticle(const Vector3 &p, const Vector3 &p0, const Real mass)
{
    m_x.push_back(p);
    m_xInitial.push_back(p);
    m_x0.push_back(p0);
    m_v.push_back(Vector3::Zero());
    m_m.push_back(mass);
    m_w.push_back(mass == 0 ? 0 : 1 / mass);
}

void SimObject::addEnergy(std::shared_ptr<EnergyT> energy)
{
    m_energies.push_back(energy);
}

void SimObject::removeEnergy(std::shared_ptr<EnergyT> energy)
{
    m_energies.erase(
        std::remove_if(
            m_energies.begin(), m_energies.end(), [&energy](auto &currentEnergy) { return energy == currentEnergy; }),
        m_energies.end());
}

void SimObject::reset()
{
    m_x = m_xInitial;
    for (auto &v : m_v) {
        v.setZero();
    }
}

void SimObject::scaleMasses(const Real scaleFactor)
{
    for (int i = 0; i < m_m.size(); ++i) {
        if (m_m[i] != 0) {
            m_m[i] *= scaleFactor;
            m_w[i] /= scaleFactor;
        }
    }
}

void RigidBodyGroup::addBody(const Vector3 &p,
                             const Vector3 &p0,
                             const Real mass,
                             const Vector3 &theta,
                             const Vector3 &theta0,
                             const Matrix3 &localInertia,
                             const TriangleMesh &mesh)
{
    m_x.push_back(p);
    m_xInitial.push_back(p);
    m_x0.push_back(p0);
    m_v.push_back(Vector3::Zero());
    m_m.push_back(mass);
    m_w.push_back(mass == 0 ? 0 : 1 / mass);
    m_theta.push_back(theta);
    m_thetaInitial.push_back(theta);
    m_omegas.push_back(Vector3::Zero());
    m_theta0.push_back(theta0);
    m_localInertia.push_back(localInertia);
    m_localInertiaInv.push_back(localInertia.inverse());
    const auto angle = theta.norm();
    m_rotationMatrix.push_back(axisAngleToRotMatrix(theta));
    m_inertia.emplace_back();
    m_inertiaInv.emplace_back();
    updateInertia(m_inertia.size() - 1);
    m_visualMesh.push_back(mesh);
}

void RigidBodyGroup::addEnergy(std::shared_ptr<EnergyT> energy)
{
    m_energies.push_back(energy);
}

void RigidBodyGroup::removeEnergy(std::shared_ptr<EnergyT> energy)
{
    m_energies.erase(
        std::remove_if(
            m_energies.begin(), m_energies.end(), [&energy](auto &currentEnergy) { return energy == currentEnergy; }),
        m_energies.end());
}

void RigidBodyGroup::reset()
{
    m_x = m_xInitial;
    for (auto &v : m_v) {
        v.setZero();
    }
    m_theta = m_thetaInitial;
    updateRotationMatrices();
    updateInertias();
    for (auto &omega : m_omegas) {
        omega.setZero();
    }
}

void RigidBodyGroup::scaleMasses(const Real scaleFactor)
{
    for (int i = 0; i < m_m.size(); ++i) {
        if (m_m[i] != 0) {
            m_m[i] *= scaleFactor;
            m_w[i] /= scaleFactor;
            m_localInertia[i] *= scaleFactor;
            m_localInertiaInv[i] /= scaleFactor;
        }
    }
}

void RigidBodyGroup::updateInertia(const int i)
{
    m_inertia[i] = m_rotationMatrix[i] * m_localInertia[i] * m_rotationMatrix[i].transpose();
    m_inertiaInv[i] = m_rotationMatrix[i] * m_localInertiaInv[i] * m_rotationMatrix[i].transpose();
}

void RigidBodyGroup::updateInertias()
{
    for (int i = 0; i < nBodies(); ++i) {
        updateInertia(i);
    }
}

void RigidBodyGroup::updateRotationMatrix(const int i)
{
    m_rotationMatrix[i] = axisAngleToRotMatrix(m_theta[i]);
}

void RigidBodyGroup::updateRotationMatrices()
{
    for (int i = 0; i < nBodies(); ++i) {
        updateRotationMatrix(i);
    }
}

void RigidBodyGroup::updateTheta(const int i)
{
    m_theta[i] = rotMatrixToAxisAngle(m_rotationMatrix[i]);
}

void RigidBodyGroup::updateThetas()
{
    for (int i = 0; i < nBodies(); ++i) {
        updateTheta(i);
    }
}

}  // namespace spg