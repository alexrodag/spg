#include <spg/sim/rigidBodyGroup.h>
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
    const Eigen::AngleAxis<spg::Real> axisAngleOrientation(R);
    return axisAngleOrientation.axis() * axisAngleOrientation.angle();
}
}  // namespace

namespace spg
{

void RigidBodyGroup::getPositions(VectorX &pos, int offsetIndex) const
{
    const int nbodies = size();
    for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
        const int startOffset = bodyIdx * 6 + offsetIndex;
        pos.segment<3>(startOffset) = m_x[bodyIdx];
        pos.segment<3>(startOffset + 3) = m_theta[bodyIdx];
    }
}

void RigidBodyGroup::getVelocities(VectorX &vel, int offsetIndex) const
{
    const int nbodies = size();
    for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
        const int startOffset = bodyIdx * 6 + offsetIndex;
        vel.segment<3>(startOffset) = m_v[bodyIdx];
        vel.segment<3>(startOffset + 3) = m_omega[bodyIdx];
    }
}

void RigidBodyGroup::setPositions(const VectorX &pos, int offsetIndex)
{
    const int nbodies = size();
    for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
        const int startOffset = bodyIdx * 6 + offsetIndex;
        m_x[bodyIdx] = pos.segment<3>(startOffset);
        m_theta[bodyIdx] = pos.segment<3>(startOffset + 3);
    }
}

void RigidBodyGroup::setVelocities(const VectorX &vel, int offsetIndex)
{
    const int nbodies = size();
    for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
        const int startOffset = bodyIdx * 6 + offsetIndex;
        m_v[bodyIdx] = vel.segment<3>(startOffset);
        m_omega[bodyIdx] = vel.segment<3>(startOffset + 3);
    }
}

void RigidBodyGroup::updateElementPositionFromDx(const Vector<6> &dx, int elementId)
{
    // Linear part
    m_x[elementId] += dx.segment<3>(0);
    // Angular part
    // Incremental compositions of the rotation
    if (const Real dxNorm = dx.segment<3>(3).norm(); dxNorm != 0) {
        m_rotationMatrix[elementId] =
            Eigen::AngleAxis<spg::Real>(dxNorm, dx.segment<3>(3) / dxNorm).toRotationMatrix() *
            m_rotationMatrix[elementId];
    }
    updateTheta(elementId);
    updateInertia(elementId);
}

void RigidBodyGroup::updatePositionsFromDx(const VectorX &dx, int offsetIndex)
{
    const int nbodies = size();
    for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
        const int startOffset = bodyIdx * 6 + offsetIndex;
        // Linear part
        m_x[bodyIdx] += dx.segment<3>(startOffset);
        // Angular part
        // Incremental compositions of the rotation
        if (const Real dxNorm = dx.segment<3>(startOffset + 3).norm(); dxNorm != 0) {
            m_rotationMatrix[bodyIdx] =
                Eigen::AngleAxis<spg::Real>(dxNorm, dx.segment<3>(startOffset + 3) / dxNorm).toRotationMatrix() *
                m_rotationMatrix[bodyIdx];
        }
    }
    updateThetas();
    updateInertias();
}

void RigidBodyGroup::integrateVelocities(Real dt)
{
    const int nbodies = size();
    for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
        // Linear part
        m_x[bodyIdx] += m_v[bodyIdx] * dt;
        // Angular part
        if (const auto omegaNorm = m_omega[bodyIdx].norm(); omegaNorm != 0) {
            // Incremental compositions of the rotation
            m_rotationMatrix[bodyIdx] =
                Eigen::AngleAxis<spg::Real>(omegaNorm * dt, m_omega[bodyIdx] / omegaNorm).toRotationMatrix() *
                m_rotationMatrix[bodyIdx];
        }
    }
    updateThetas();
    updateInertias();
}

void RigidBodyGroup::computeIntegratedVelocities(const VectorX &oldPos, int offsetIndex, Real invdt)
{
    const int nbodies = size();
    for (int bodyIdx = 0; bodyIdx < nbodies; ++bodyIdx) {
        const int startOffset = bodyIdx * 6 + offsetIndex;
        // Linear part
        m_v[bodyIdx] = (m_x[bodyIdx] - oldPos.segment<3>(startOffset)) * invdt;
        // Angular part
        Eigen::AngleAxis<spg::Real> deltaAxisAngle(
            m_rotationMatrix[bodyIdx] * axisAngleToRotMatrix(oldPos.segment<3>(startOffset + 3)).transpose());
        m_omega[bodyIdx] = deltaAxisAngle.angle() * deltaAxisAngle.axis() * invdt;
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
    m_omega.push_back(Vector3::Zero());
    m_theta0.push_back(theta0);
    m_localInertia.push_back(localInertia);
    m_localInertiaInv.push_back(localInertia.inverse());
    const auto angle = theta.norm();
    m_rotationMatrix.push_back(axisAngleToRotMatrix(theta));
    m_inertia.emplace_back();
    m_inertiaInv.emplace_back();
    updateInertia(static_cast<int>(m_inertia.size()) - 1);
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
    for (auto &omega : m_omega) {
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
    for (int i = 0; i < size(); ++i) {
        updateInertia(i);
    }
}

void RigidBodyGroup::updateRotationMatrix(const int i)
{
    m_rotationMatrix[i] = axisAngleToRotMatrix(m_theta[i]);
}

void RigidBodyGroup::updateRotationMatrices()
{
    for (int i = 0; i < size(); ++i) {
        updateRotationMatrix(i);
    }
}

void RigidBodyGroup::updateTheta(const int i)
{
    m_theta[i] = rotMatrixToAxisAngle(m_rotationMatrix[i]);
}

void RigidBodyGroup::updateThetas()
{
    for (int i = 0; i < size(); ++i) {
        updateTheta(i);
    }
}

}  // namespace spg