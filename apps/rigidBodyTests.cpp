#include <iostream>
#include <vector>
#include <polyscope/polyscope.h>
#include <polyscope/combining_hash_functions.h>
#include <polyscope/messages.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/volume_mesh.h>
#include <polyscope/pick.h>
#include <Eigen/Dense>
#include <numbers>

#include <spg/types.h>
#include <spg/geom/triangleMesh.h>
#include <spg/geom/io/obj.h>

/////////////////////////////////////////
// Polyscope and custom utils for picking
/////////////////////////////////////////

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

///////////////
// Main program
///////////////

spg::Matrix3 rodrigues(const spg::Vector3 &v)
{
    const auto vNorm = v.norm();
    if (vNorm == 0) {
        return spg::Matrix3::Identity();
    }
    return Eigen::AngleAxis<spg::Real>(vNorm, v / vNorm).toRotationMatrix();
}
spg::Matrix3 skew(const spg::Vector3 &v)
{
    spg::Matrix3 vSkew;
    vSkew << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
    return vSkew;
}

spg::Matrix3 dRdvGallego(const spg::Matrix3 &R, const spg::Vector3 &v)
{
    std::array<spg::Matrix3, 3> dRdviTimesRt;
    const spg::Real vSqNormInv = 1 / v.squaredNorm();
    const spg::Matrix3 I = spg::Matrix3::Identity();
    for (int i = 0; i < 3; ++i) {
        dRdviTimesRt[i] = vSqNormInv * (v[i] * skew(v) + skew(v.cross((I - R) * I.col(i))));
        // It should be multiplied by R also, but since we are then multipying by R^T, they cancel out
        // std::cout << "dRdvGallego" << i << "\n" << dRdviTimesRt[i] << "\n\n";
    }
    return spg::Matrix3{{dRdviTimesRt[0](2, 1), dRdviTimesRt[1](2, 1), dRdviTimesRt[2](2, 1)},
                        {dRdviTimesRt[0](0, 2), dRdviTimesRt[1](0, 2), dRdviTimesRt[2](0, 2)},
                        {dRdviTimesRt[0](1, 0), dRdviTimesRt[1](1, 0), dRdviTimesRt[2](1, 0)}};
}

struct RigidBody {
    enum class Mode { IncrementalRodrigues, IncrementalSmallAngleRodrigues, Gallego, ETHZ, Quaternion };
    spg::Vector3 pos;
    spg::Vector3 orientation;
    Eigen::Quaternion<spg::Real> Q;
    spg::Vector3 vel;
    spg::Vector3 omega;
    spg::Matrix3 localInertiaTensor;
    spg::Matrix3 localInertiaTensorInv;
    spg::Matrix3 inertiaTensor;
    spg::Matrix3 inertiaTensorInv;
    spg::Matrix3 rotationMatrix;
    spg::TriangleMesh visualMesh;
    Mode mode{Mode::IncrementalRodrigues};

    RigidBody(const spg::TriangleMesh &mesh)
        : visualMesh(mesh)
    {
    }

    void updateInertiaTensor()
    {
        inertiaTensor = rotationMatrix * localInertiaTensor * rotationMatrix.transpose();
        inertiaTensorInv = rotationMatrix * localInertiaTensorInv * rotationMatrix.transpose();
    }

    void integrateVelocity(spg::Real dt)
    {
        vel = vel;
        omega += dt * inertiaTensorInv * (-omega.cross(inertiaTensor * omega));
    }

    void integratePositions(spg::Real dt)
    {
        constexpr spg::Real pi = 3.14159265358979323846;
        pos = pos + dt * vel;
        if (mode == Mode::IncrementalRodrigues) {
            if (const auto omegaNorm = omega.norm(); omegaNorm != 0) {
                // Incremental compositions of the rotation
                Eigen::AngleAxis<spg::Real> axisAngleOrientation(
                    Eigen::AngleAxis<spg::Real>(omegaNorm * dt, omega / omegaNorm).toRotationMatrix() * rotationMatrix);
                orientation = axisAngleOrientation.angle() * axisAngleOrientation.axis();
                rotationMatrix = axisAngleOrientation.toRotationMatrix();
            }
        }
        if (mode == Mode::IncrementalSmallAngleRodrigues) {
            Eigen::AngleAxis<spg::Real> axisAngleOrientation((spg::Matrix3::Identity() + skew(omega * dt)) *
                                                             rotationMatrix);
            orientation = axisAngleOrientation.angle() * axisAngleOrientation.axis();
            rotationMatrix = axisAngleOrientation.toRotationMatrix();
        }
        if (mode == Mode::ETHZ) {
            const auto skewTheta = skew(orientation);
            const auto angle = orientation.norm();
            const spg::Matrix3 J =
                spg::Matrix3::Identity() - 0.5 * skewTheta +
                skewTheta * skewTheta * (1 - 0.5 * angle * sin(angle) / (1 - cos(angle))) / (angle * angle);
            orientation += J * (omega * dt);
            auto newAngle = orientation.norm();
            const spg::Vector3 newAxis = orientation / newAngle;
            while (newAngle > pi) {
                newAngle -= 2 * pi;
            }
            while (newAngle < -pi) {
                newAngle += 2 * pi;
            }
            orientation = newAxis * newAngle;
            rotationMatrix = Eigen::AngleAxis<spg::Real>(newAngle, newAxis).toRotationMatrix();
        }
        if (mode == Mode::Gallego) {
            const spg::Matrix3 J = dRdvGallego(rodrigues(orientation), orientation).inverse();
            orientation += J * (omega * dt);
            auto newAngle = orientation.norm();
            const spg::Vector3 newAxis = orientation / newAngle;
            while (newAngle > pi) {
                newAngle -= 2 * pi;
            }
            while (newAngle < -pi) {
                newAngle += 2 * pi;
            }
            orientation = newAxis * newAngle;
            rotationMatrix = Eigen::AngleAxis<spg::Real>(newAngle, newAxis).toRotationMatrix();
        }
        if (mode == Mode::Quaternion) {
            // Last and first rows are swapped since Eigen Quaternion coeffs have w as last one
            Eigen::Matrix<spg::Real, 4, 3> J{{Q.w(), Q.z(), -Q.y()},  //
                                             {-Q.z(), Q.w(), Q.x()},
                                             {Q.y(), -Q.x(), Q.w()},
                                             {-Q.x(), -Q.y(), -Q.z()}};
            J *= 0.5;
            Q.coeffs() += J * (omega * dt);
            Q.normalize();
            rotationMatrix = Q.toRotationMatrix();
        }
    }

    void step(spg::Real dt)
    {
        integrateVelocity(dt);
        integratePositions(dt);
        updateInertiaTensor();
    }

    void prepare()
    {
        rotationMatrix =
            orientation.norm() != 0
                ? Eigen::AngleAxis<spg::Real>(orientation.norm(), orientation.normalized()).toRotationMatrix()
                : spg::Matrix3::Identity();
        updateInertiaTensor();
        const auto vNorm = orientation.norm();
        Q = Eigen::Quaternion<spg::Real>(Eigen::AngleAxis<spg::Real>(vNorm, orientation / vNorm));
    }
};
std::tuple<std::vector<spg::Vector3>, std::vector<spg::Int3>> unitCubeGeometry()
{
    std::vector<spg::Vector3> vertices{{-0.5, -0.5, -0.5},
                                       {0.5, -0.5, -0.5},
                                       {-0.5, 0.5, -0.5},
                                       {0.5, 0.5, -0.5},
                                       {-0.5, -0.5, 0.5},
                                       {0.5, -0.5, 0.5},
                                       {-0.5, 0.5, 0.5},
                                       {0.5, 0.5, 0.5}};
    std::vector<spg::Int3> faces{{0, 2, 3},
                                 {0, 3, 1},
                                 {0, 4, 6},
                                 {0, 6, 2},
                                 {0, 1, 5},
                                 {0, 5, 4},
                                 {7, 6, 4},
                                 {7, 4, 5},
                                 {7, 3, 2},
                                 {7, 2, 6},
                                 {7, 5, 1},
                                 {7, 1, 3}};
    return {vertices, faces};
}
int main()
{
    try {
        int currentPickedParticleId = -1;
        float currentPickedDepth = -1;
        spg::Vector3 currentPickedPos;
        spg::Vector3 pickAnchorPos;

        // Options
        // polyscope::view::windowWidth = 600;
        // polyscope::options::maxFPS = -1;
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
        // polyscope::view::setFrontDir(polyscope::FrontDir::NegZFront);
        // polyscope::view::setNavigateStyle(polyscope::NavigateStyle::Free);

        // Initialize polyscope
        polyscope::init();

        // Runtime configurable values
        bool run{false};
        float dt{0.01f};
        int substeps{1};

        auto [vertices, faces] = unitCubeGeometry();
        spg::Real width{1}, height{3}, depth{0.2}, mass{1};
        spg::Matrix3 transform;
        transform.setZero();
        transform(0, 0) = width;
        transform(1, 1) = height;
        transform(2, 2) = depth;
        for (auto &v : vertices) {
            v = transform * v;
        }
        spg::TriangleMesh mesh(vertices, vertices, faces);
        std::vector<RigidBody> rbs;
        rbs.push_back(RigidBody(mesh));
        auto &rigidBody = rbs.front();

        rbs.front().localInertiaTensor.setZero();
        rbs.front().localInertiaTensor(0, 0) = 1.0 / 12.0 * mass * (height * height + depth * depth);
        rbs.front().localInertiaTensor(1, 1) = 1.0 / 12.0 * mass * (width * width + depth * depth);
        rbs.front().localInertiaTensor(2, 2) = 1.0 / 12.0 * mass * (width * width + height * height);
        rbs.front().localInertiaTensorInv = rbs.front().localInertiaTensor.inverse();

        rbs.front().pos = {0, 0, 0};
        rbs.front().vel = {0, 0, 0};
        rbs.front().orientation = {0, 0, 0.001};
        rbs.front().orientation = {0.2, 0.3, 0.401};
        rbs.front().omega = {1, 0.000001, 0.000001};
        rbs.front().omega = {1, 0.7, 0.3};
        /* rbs.front().omega = {0.000001, 1, 0.000001}; */
        /* rbs.front().omega = {0.000001, 0.000001, 1}; */
        /* rbs.front().omega = {0, 0, 0}; */
        /* rbs.front().omega = {0.3, 0.2, 0.1}; */
        rbs.push_back(rbs.front());
        rbs.back().mode = RigidBody::Mode::ETHZ;
        rbs.push_back(rbs.front());
        rbs.back().mode = RigidBody::Mode::Gallego;
        rbs.push_back(rbs.front());
        rbs.back().mode = RigidBody::Mode::Quaternion;
        rbs.push_back(rbs.front());
        rbs.back().mode = RigidBody::Mode::IncrementalSmallAngleRodrigues;
        std::vector<polyscope::SurfaceMesh *> rigidBodySM;
        int id = 0;
        for (auto &rb : rbs) {
            rb.prepare();
            rigidBodySM.push_back(polyscope::registerSurfaceMesh(
                "RB" + std::to_string(id), rb.visualMesh.vertices(), rb.visualMesh.faces()));
            ++id;
        }
        auto l_callback = [&]() {
            ImGui::PushItemWidth(200);
            ImGui::Checkbox("Simulate [Space]", &run);
            ImGui::Dummy(ImVec2(0.0f, 20.0f));
            if (ImGui::InputFloat("dt", &dt, 0, 0, "%.6f")) {
                dt = std::max(dt, 1e-6f);
            }
            if (ImGui::InputInt("substeps", &substeps)) {
                substeps = std::max(substeps, 1);
            }
            if (ImGui::Button("Reset solvers [r]")) {
            }
            auto l_rbStep = [&]() {
                for (int i = 0; i < substeps; ++i) {
                    for (auto &rb : rbs) {
                        rb.step(dt / substeps);
                    }
                }
                for (int i = 0; i < rbs.size(); ++i) {
                    auto &rb = rbs[i];
                    auto &rbmat = rb.rotationMatrix;
                    auto &pos = rb.pos;
                    glm::mat4 mat(rbmat(0, 0),
                                  rbmat(0, 1),
                                  rbmat(0, 2),
                                  pos(0),
                                  rbmat(1, 0),
                                  rbmat(1, 1),
                                  rbmat(1, 2),
                                  pos(1),
                                  rbmat(2, 0),
                                  rbmat(2, 1),
                                  rbmat(2, 2),
                                  pos(2),
                                  0,
                                  0,
                                  0,
                                  1);
                    mat = glm::transpose(mat);
                    rigidBodySM[i]->setTransform(mat);
                    // std::cout << "O"+std::to_string(i)+" " << rigidBody.orientation.transpose() << "\n";
                    // std::cout << "w"+std::to_string(i)+" " << rigidBody.omega.transpose() << "\n";
                    // std::cout << "L" + std::to_string(i) + " " << (rb.inertiaTensor * rb.omega).transpose() << "\n";
                }
            };

            if (ImGui::Button("Single Step [s]")) {
                if (!run) {
                    l_rbStep();
                }
            }
            if (run) {
                l_rbStep();
            }

            if (ImGui::Button("Center camera")) {
                polyscope::view::resetCameraToHomeView();
            }

            ImGui::Dummy(ImVec2(0.0f, 20.0f));
            ImGui::Text("Mid-click and drag on a vertex to pull it");
            ImGui::Dummy(ImVec2(0.0f, 20.0f));

            ImGui::PopItemWidth();

            if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_Space)) {
                run = !run;
            }

            if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_S)) {
                if (!run) {
                    l_rbStep();
                }
            }
        };
        polyscope::state::userCallback = l_callback;

        // Show the gui
        polyscope::show();
    } catch (const std::exception e) {
        std::cerr << "Exception thrown: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}