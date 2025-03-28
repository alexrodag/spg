#include <spg/sim/simObject/particleGroup.h>
#include <spg/sim/energy/energy.h>

namespace spg
{

void ParticleGroup::getPositions(VectorX &pos, int offsetIndex) const
{
    const int nparticles = size();
    for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
        const int startOffset = particleIdx * 3 + offsetIndex;
        pos.segment<3>(startOffset) = m_x[particleIdx];
    }
}

void ParticleGroup::getVelocities(VectorX &vel, int offsetIndex) const
{
    const int nparticles = size();
    for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
        const int startOffset = particleIdx * 3 + offsetIndex;
        vel.segment<3>(startOffset) = m_v[particleIdx];
    }
}

void ParticleGroup::setPositions(const VectorX &pos, int offsetIndex)
{
    const int nparticles = size();
    for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
        const int startOffset = particleIdx * 3 + offsetIndex;
        m_x[particleIdx] = pos.segment<3>(startOffset);
    }
}

void ParticleGroup::setVelocities(const VectorX &vel, int offsetIndex)
{
    const int nparticles = size();
    for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
        const int startOffset = particleIdx * 3 + offsetIndex;
        m_v[particleIdx] = vel.segment<3>(startOffset);
    }
}

void ParticleGroup::updateElementPositionFromDx(const Vector3 &dx, int elementId)
{
    m_x[elementId] += dx;
}

void ParticleGroup::updatePositionsFromDx(const VectorX &dx, int offsetIndex)
{
    const int nparticles = size();
    for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
        const int startOffset = particleIdx * 3 + offsetIndex;
        m_x[particleIdx] += dx.segment<3>(startOffset);
    }
}

void ParticleGroup::integrateVelocities(Real dt)
{
    const int nparticles = size();
    for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
        m_x[particleIdx] += m_v[particleIdx] * dt;
    }
}

void ParticleGroup::computeIntegratedVelocities(const VectorX &oldPos, int offsetIndex, Real invdt)
{
    const int nparticles = size();
    for (int particleIdx = 0; particleIdx < nparticles; ++particleIdx) {
        const int startOffset = particleIdx * 3 + offsetIndex;
        m_v[particleIdx] = (m_x[particleIdx] - oldPos.segment<3>(startOffset)) * invdt;
    }
}

void ParticleGroup::addParticle(const Vector3 &p, const Vector3 &p0, const Real mass)
{
    m_x.push_back(p);
    m_xInitial.push_back(p);
    m_x0.push_back(p0);
    m_v.push_back(Vector3::Zero());
    m_m.push_back(mass);
    m_w.push_back(mass == 0 ? 0 : 1 / mass);
}

void ParticleGroup::addEnergy(std::shared_ptr<EnergyT> energy)
{
    m_energies.push_back(energy);
}

void ParticleGroup::removeEnergy(std::shared_ptr<EnergyT> energy)
{
    m_energies.erase(
        std::remove_if(
            m_energies.begin(), m_energies.end(), [&energy](auto &currentEnergy) { return energy == currentEnergy; }),
        m_energies.end());
}

void ParticleGroup::reset()
{
    m_x = m_xInitial;
    for (auto &v : m_v) {
        v.setZero();
    }
}

void ParticleGroup::scaleMasses(const Real scaleFactor)
{
    for (int i = 0; i < static_cast<int>(m_m.size()); ++i) {
        if (m_m[i] != 0) {
            m_m[i] *= scaleFactor;
            m_w[i] /= scaleFactor;
        }
    }
}

}  // namespace spg