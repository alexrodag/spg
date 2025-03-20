#include <iostream>
#include <vector>
#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include <polyscope/pick.h>

#include <spg/types.h>
#include <spg/sim/simObject/particleGroup.h>
#include <spg/sim/energy/energy.h>
#include <spg/sim/energy/springEnergy.h>
#include <spg/sim/energy/springContinuumEnergy.h>
#include <spg/sim/energy/discreteBendingEnergy.h>
#include <spg/sim/energy/baraffWitkinBendingEnergy.h>
#include <spg/sim/energy/quadraticBendingEnergy.h>
#include <spg/sim/energy/membraneStvkEnergy.h>
#include <spg/sim/energy/membraneBaraffWitkinEnergy.h>
#include <spg/sim/energy/membraneChoiEnergy.h>
#include <spg/sim/energy/stableNeoHookeanEnergy.h>

spg::ParticleGroup createStvkObject()
{
    spg::ParticleGroup obj;
    obj.addParticle({0, 0, 0}, {0, 0, 0}, 1);
    obj.addParticle({1, 0, 0}, {1, 0, 0}, 1);
    obj.addParticle({0, 1, 0}, {0, 1, 0}, 1);
    auto stvkEnergy = std::make_shared<spg::MembraneStvkEnergy>();
    stvkEnergy->addStencil({0, 1, 2}, 1, 0);
    obj.addEnergy(stvkEnergy);
    stvkEnergy->preparePrecomputations(obj);
    return obj;
}

spg::ParticleGroup createBaraffWitkinObject()
{
    spg::ParticleGroup obj;
    obj.addParticle({0, 0, 0}, {0, 0, 0}, 1);
    obj.addParticle({1, 0, 0}, {1, 0, 0}, 1);
    obj.addParticle({0, 1, 0}, {0, 1, 0}, 1);
    auto stvkEnergy = std::make_shared<spg::MembraneBaraffWitkinEnergy>();
    stvkEnergy->addStencil({0, 1, 2}, 1, 1, 1);
    obj.addEnergy(stvkEnergy);
    stvkEnergy->preparePrecomputations(obj);
    return obj;
}

spg::ParticleGroup createChoiObject()
{
    spg::ParticleGroup obj;
    obj.addParticle({0, 0, 0}, {0, 0, 0}, 1);
    obj.addParticle({1, 0, 0}, {1, 0, 0}, 1);
    obj.addParticle({0, 1, 0}, {0, 1, 0}, 1);
    auto stvkEnergy = std::make_shared<spg::MembraneChoiEnergy>();
    stvkEnergy->addStencil({0, 1, 2}, 1, 1, 1);
    obj.addEnergy(stvkEnergy);
    stvkEnergy->preparePrecomputations(obj);
    return obj;
}

spg::ParticleGroup createDiscreteBendingObject()
{
    spg::ParticleGroup obj;
    obj.addParticle({1, 0, 1}, {1, -1, 0}, 1);
    obj.addParticle({1, 0, -1}, {1, 1, 0}, 1);
    obj.addParticle({0, 0, 0}, {0, 0, 0}, 1);
    obj.addParticle({2, 0, 0}, {2, 0, 0}, 1);
    auto bendingEnergy = std::make_shared<spg::DiscreteBendingEnergy>();
    bendingEnergy->addStencil({0, 1, 2, 3}, 0, 1);
    obj.addEnergy(bendingEnergy);
    bendingEnergy->preparePrecomputations(obj);
    return obj;
}

spg::ParticleGroup createBWBendingObject()
{
    spg::ParticleGroup obj;
    obj.addParticle({1, 0, 1}, {1, -1, 0}, 1);
    obj.addParticle({1, 0, -1}, {1, 1, 0}, 1);
    obj.addParticle({0, 0, 0}, {0, 0, 0}, 1);
    obj.addParticle({2, 0, 0}, {2, 0, 0}, 1);
    auto bendingEnergy = std::make_shared<spg::BaraffWitkinBendingEnergy>();
    bendingEnergy->addStencil({0, 1, 2, 3}, 0, 1);
    obj.addEnergy(bendingEnergy);
    bendingEnergy->preparePrecomputations(obj);
    return obj;
}

spg::ParticleGroup createQuadraticBendingObject()
{
    spg::ParticleGroup obj;
    obj.addParticle({1, 0, 1}, {1, -1, 0}, 1);
    obj.addParticle({1, 0, -1}, {1, 1, 0}, 1);
    obj.addParticle({0, 0, 0}, {0, 0, 0}, 1);
    obj.addParticle({2, 0, 0}, {2, 0, 0}, 1);
    auto bendingEnergy = std::make_shared<spg::QuadraticBendingEnergy>();
    bendingEnergy->addStencil({0, 1, 2, 3}, 0, 1);
    obj.addEnergy(bendingEnergy);
    bendingEnergy->preparePrecomputations(obj);
    return obj;
}

spg::ParticleGroup createSpringObject()
{
    spg::ParticleGroup obj;
    obj.addParticle({0, 0, 0}, {0, 0, 0}, 1);
    obj.addParticle({1, 0, 0}, {1, 0, 0}, 1);
    auto springEnergy = std::make_shared<spg::SpringEnergy>();
    springEnergy->addStencil({0, 1}, 1, 1);
    obj.addEnergy(springEnergy);
    springEnergy->preparePrecomputations(obj);
    return obj;
}

spg::ParticleGroup createSpringContinuumObject()
{
    spg::ParticleGroup obj;
    obj.addParticle({0, 0, 0}, {0, 0, 0}, 1);
    obj.addParticle({1, 0, 0}, {1, 0, 0}, 1);
    auto springContinuumEnergy = std::make_shared<spg::SpringContinuumEnergy>();
    springContinuumEnergy->addStencil({0, 1}, 1, 1);
    obj.addEnergy(springContinuumEnergy);
    springContinuumEnergy->preparePrecomputations(obj);
    return obj;
}

spg::ParticleGroup createStableNeoHookeanObject()
{
    spg::ParticleGroup obj;
    obj.addParticle({0, 0, 0}, {0, 0, 0}, 1);
    obj.addParticle({1, 0, 0}, {1, 0, 0}, 1);
    obj.addParticle({0, 1, 0}, {0, 1, 0}, 1);
    obj.addParticle({0, 0, 1}, {0, 0, 1}, 1);
    auto stableNeoHookeanEnergy = std::make_shared<spg::StableNeoHookeanEnergy>();
    stableNeoHookeanEnergy->addStencil({0, 1, 2, 3}, 10, 0.2);
    obj.addEnergy(stableNeoHookeanEnergy);
    stableNeoHookeanEnergy->preparePrecomputations(obj);
    return obj;
}

// Polyscope utils

glm::vec3 screenCoordsToWorldPosition(glm::vec2 screenCoords, float &depth)
{
    glm::mat4 view = polyscope::view::getCameraViewMatrix();
    glm::mat4 viewInv = glm::inverse(view);
    glm::mat4 proj = polyscope::view::getCameraPerspectiveMatrix();
    glm::mat4 projInv = glm::inverse(proj);
    if (depth < 0.) {
        int xInd, yInd;
        std::tie(xInd, yInd) = polyscope::view::screenCoordsToBufferInds(screenCoords);
        // query the depth buffer to get depth
        polyscope::render::FrameBuffer *sceneFramebuffer = polyscope::render::engine->sceneBuffer.get();
        depth = sceneFramebuffer->readDepth(xInd, polyscope::view::bufferHeight - yInd);
        if (depth == 1.) {
            // if we didn't hit anything in the depth buffer, just return infinity
            float inf = std::numeric_limits<float>::infinity();
            return glm::vec3{inf, inf, inf};
        }
    }

    // convert depth to world units
    glm::vec2 screenPos{screenCoords.x / static_cast<float>(polyscope::view::windowWidth),
                        1.f - screenCoords.y / static_cast<float>(polyscope::view::windowHeight)};
    float z = depth * 2.0f - 1.0f;
    glm::vec4 clipPos = glm::vec4(screenPos * 2.0f - 1.0f, z, 1.0f);
    glm::vec4 viewPos = projInv * clipPos;
    viewPos /= viewPos.w;

    glm::vec4 worldPos = viewInv * viewPos;
    worldPos /= worldPos.w;

    return glm::vec3(worldPos);
}

int main(int argc, char **argv)
{
    // Stuff for picking
    int currentPickedParticleId = -1;
    float currentPickedDepth = -1;
    spg::Vector3 currentPickedPos;

    // Options
    // polyscope::view::windowWidth = 600;
    // polyscope::options::maxFPS = -1;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    // polyscope::view::setFrontDir(polyscope::FrontDir::NegZFront);
    //  polyscope::view::setNavigateStyle(polyscope::NavigateStyle::Free);
    //   Initialize polyscope
    polyscope::init();

    auto obj = createStvkObject();
    int objType = 0;
    float regularizerScale = 1e-5;

    // Register obj
    auto *pc = polyscope::registerPointCloud("simObj", obj.positions());
    pc->setPointRadius(0.1, false);
    auto l_callback = [&]() {
        auto l_updateView = [&obj](const spg::MatrixX &forces = {}, const spg::MatrixX &hessianWeightedForces = {}) {
            auto *pc = polyscope::getPointCloud("simObj");
            pc->updatePointPositions(obj.positions());
            if (forces.size() != 0) {
                auto *vc = pc->addVectorQuantity("Forces", forces, polyscope::VectorType::AMBIENT);
                vc->setVectorRadius(0.1, false);
                // vc->setEnabled(true);
            } else {
                pc->removeQuantity("Forces");
            }
            if (hessianWeightedForces.size() != 0) {
                auto *vc = pc->addVectorQuantity(
                    "Hessian weighted forces", hessianWeightedForces, polyscope::VectorType::AMBIENT);
                vc->setVectorRadius(0.1, false);
                // vc->setEnabled(true);
            } else {
                pc->removeQuantity("Hessian weighted forces");
            }
        };

        auto l_getForces = [&obj, regularizerScale]() {
            const spg::VectorX grad = obj.energies().front()->energyGradientGeneric(0, obj);
            spg::MatrixX hess = obj.energies().front()->energyHessianGeneric(0, obj);
            spg::MatrixX forces;
            forces.resize(obj.size(), 3);
            for (int i = 0; i < forces.rows(); ++i) {
                forces.row(i) = -grad.segment<3>(i * 3);
            }
            for (int i = 0; i < obj.size(); ++i) {
                hess.block<3, 3>(i * 3, i * 3) += spg::Matrix3::Identity() * regularizerScale;
            }
            spg::VectorX hessianInverseTimesGrad = hess.inverse() * grad;
            spg::MatrixX hessianWeighedForces;
            hessianWeighedForces.resize(obj.size(), 3);
            for (int i = 0; i < hessianWeighedForces.rows(); ++i) {
                hessianWeighedForces.row(i) = -hessianInverseTimesGrad.segment<3>(i * 3);
            }
            return std::tuple<spg::MatrixX, spg::MatrixX>(forces, hessianWeighedForces);
        };
        auto l_registerObj = [&obj]() {
            polyscope::registerPointCloud("simObj", obj.positions())->setPointRadius(0.1, false);
        };
        auto l_unregisterObj = []() { polyscope::removePointCloud("simObj"); };
        ImGui::PushItemWidth(200);
        const char *objTypeNames[] = {
            "StVK",
            "BaraffWitkin membrane",
            "Choi membrane",
            "Bending",
            "BWBending",
            "QuadraticBending",
            "Spring",
            "ContinuumSpring",
            "StableNeoHookean",
        };

        if (ImGui::ListBox(std::string("Energy type").c_str(), &objType, objTypeNames, IM_ARRAYSIZE(objTypeNames))) {
            l_unregisterObj();
            if (objType == 0) {
                obj = createStvkObject();
            } else if (objType == 1) {
                obj = createBaraffWitkinObject();
            } else if (objType == 2) {
                obj = createChoiObject();
            } else if (objType == 3) {
                obj = createDiscreteBendingObject();
            } else if (objType == 4) {
                obj = createBWBendingObject();
            } else if (objType == 5) {
                obj = createQuadraticBendingObject();
            } else if (objType == 6) {
                obj = createSpringObject();
            } else if (objType == 7) {
                obj = createSpringContinuumObject();
            } else if (objType == 8) {
                obj = createStableNeoHookeanObject();
            }
            l_registerObj();
            const auto [forces, hessianWeightedForces] = l_getForces();
            l_updateView(forces, hessianWeightedForces);
        }

        if (ImGui::Button("Reset")) {
            obj.reset();
            const auto [forces, hessianWeightedForces] = l_getForces();
            l_updateView(forces, hessianWeightedForces);
        }
        if (ImGui::SliderFloat("Regularizer scale", &regularizerScale, 0, 10, "%.6f")) {
            const auto [forces, hessianWeightedForces] = l_getForces();
            l_updateView(forces, hessianWeightedForces);
        }

        ImGui::PopItemWidth();

        ImGuiIO &io = ImGui::GetIO();
        if (io.MouseClicked[2]) {  // if the left mouse button was clicked
            // gather values
            glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
            std::pair<polyscope::Structure *, size_t> pickPair =
                polyscope::pick::evaluatePickQuery(screenCoords.x, screenCoords.y);
            if (pickPair.first != nullptr) {
                auto pos = screenCoordsToWorldPosition(screenCoords, currentPickedDepth);
                // Not sure why it can return 1 when picking something, but it is the case, ignore the hit
                if (currentPickedDepth != 1.) {
                    currentPickedParticleId = pickPair.second;
                    currentPickedPos = {pos.x, pos.y, pos.z};
                }
            }
        } else if (io.MouseReleased[2]) {
            if (currentPickedParticleId != -1) {
                currentPickedDepth = -1;
                currentPickedParticleId = -1;
            }
        } else if (io.MouseDown[2] && (io.MouseDelta.x != 0 || io.MouseDelta.y != 0) && currentPickedDepth > 0) {
            glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
            glm::vec3 worldPos = screenCoordsToWorldPosition(screenCoords, currentPickedDepth);

            if (currentPickedParticleId != -1) {
                spg::Vector3 disp = spg::Vector3(worldPos.x, worldPos.y, worldPos.z) - currentPickedPos;
                currentPickedPos = {worldPos.x, worldPos.y, worldPos.z};
                obj.positions()[currentPickedParticleId] += disp;
                const auto [forces, hessianWeightedForces] = l_getForces();
                l_updateView(forces, hessianWeightedForces);
            }
        }
    };
    polyscope::state::userCallback = l_callback;

    // Show the gui
    polyscope::show();

    return 0;
}