#include <iostream>
#include <vector>

#include <spg/types.h>
#include <spg/sim/simObject/particleGroup.h>
#include <spg/sim/energy/springEnergy.h>
#include <spg/sim/energy/springAnchorEnergy.h>
#include <spg/sim/solver/implicitEulerBaraffWitkin.h>

int main()
{
    spg::ParticleGroup obj;
    obj.addParticle({0, 0, 0}, {0, 0, 0}, 1);
    obj.addParticle({1, 0, 0}, {1, 0, 0}, 1);
    auto springAnchorEnergy = std::make_shared<spg::SpringAnchorEnergy>();
    auto springEnergy = std::make_shared<spg::SpringEnergy>();
    springAnchorEnergy->addStencil(std::array<int, 1>{0}, obj.positions0()[0], 1e7);  // anchor the first particle
    springEnergy->addStencil(std::array<int, 2>{0, 1}, 1, 10);  // add spring between the two particles
    obj.addEnergy(springAnchorEnergy);
    obj.addEnergy(springEnergy);
    spg::solver::ImplicitEulerBaraffWitkin solver;
    solver.addObject(obj);
    const float dt = 0.1f;
    solver.setDt(dt);
    float time = 0;
    for (int i = 0; i < 100; ++i) {
        solver.step();
        time += dt;
        std::cout << "Anchored particle position at " << time
                  << "s: " << solver.particleGroups().front().positions()[0].transpose() << "\n";
        std::cout << "Free particle position at " << time
                  << "s: " << solver.particleGroups().front().positions()[1].transpose() << "\n\n";
    }
    return 0;
}