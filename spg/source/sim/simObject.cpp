#include <spg/sim/simObject.h>
#include <spg/sim/energy/energy.h>

namespace spg
{
void SimObject::addParticle(const Vector3 &p, const Vector3 &p0, Real mass)
{
    m_x.push_back(p);
    m_xInitial.push_back(p);
    m_x0.push_back(p0);
    m_v.push_back(Vector3::Zero());
    m_m.push_back(mass);
    m_w.push_back(mass == 0 ? 0 : 1 / mass);
}

void SimObject::addEnergy(std::shared_ptr<Energy> energy)
{
    m_energies.push_back(energy);
}

void SimObject::removeEnergy(std::shared_ptr<Energy> energy)
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
}  // namespace spg