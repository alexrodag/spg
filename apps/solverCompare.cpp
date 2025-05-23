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

#include <spg/types.h>
#include <spg/sim/simObject/particleGroup.h>
#include <spg/sim/simObject/rigidBodyGroup.h>
#include <spg/sim/energy/springEnergy.h>
#include <spg/sim/energy/springContinuumEnergy.h>
#include <spg/sim/energy/springSquaredContinuumEnergy.h>
#include <spg/sim/energy/springAnchorEnergy.h>
#include <spg/sim/energy/discreteBendingEnergy.h>
#include <spg/sim/energy/baraffWitkinBendingEnergy.h>
#include <spg/sim/energy/quadraticBendingEnergy.h>
#include <spg/sim/energy/membraneStvkEnergy.h>
#include <spg/sim/energy/membraneBaraffWitkinEnergy.h>
#include <spg/sim/energy/membraneChoiEnergy.h>
#include <spg/sim/energy/stableNeoHookeanEnergy.h>
#include <spg/sim/energy/stableSquaredNeoHookeanEnergy.h>
#include <spg/sim/energy/stvkEnergy.h>
#include <spg/sim/energy/rigidBodySpringAnchorEnergy.h>
#include <spg/sim/energy/rigidBodySpringEnergy.h>
#include <spg/sim/energy/rigidBodySpringSquaredEnergy.h>
#include <spg/sim/energy/rigidBodyOrientationAlignmentEnergy.h>
#include <spg/sim/solver/xpbd.h>
#include <spg/sim/solver/simplecticEuler.h>
#include <spg/sim/solver/bdf2.h>
#include <spg/sim/solver/implicitEulerBaraffWitkin.h>
#include <spg/sim/solver/implicitEulerNewtonDv.h>
#include <spg/sim/solver/implicitEulerNewtonDx.h>
#include <spg/sim/solver/implicitEulerNewtonRobust.h>
#include <spg/sim/solver/quasiStaticNewton.h>
#include <spg/sim/solver/quasiStaticNewtonRobust.h>
#include <spg/sim/solver/staticNewton.h>
#include <spg/sim/solver/vbd.h>
#include <spg/geom/triangleMesh.h>
#include <spg/geom/tetrahedralMesh.h>

///////////////////////////////////////////////////////////
// Simulation object creation methods and types for the GUI
///////////////////////////////////////////////////////////

enum class MembraneType { Stvk, BaraffWitkin, Choi };
enum class BendingType { Discrete, Quadratic, BaraffWitkin };
enum class SpringType { Normal, Continuum, SquaredContinuum };
enum class FemType { Stvk, NeoHookean, SquaredNeoHookean };

spg::ParticleGroup createRope(const float mass, const float ropeLength, const int nParticles, SpringType springType)
{
    spg::ParticleGroup pGroup;
    for (int i = 0; i < nParticles; ++i) {
        pGroup.addParticle({-ropeLength * 0.5f + i * ropeLength / (nParticles - 1),
                            0,
                            0 + (i == 0 ? spg::Real(ropeLength * 0.01) / (nParticles - 1) : 0)},
                           {-ropeLength * 0.5f + i * ropeLength / (nParticles - 1),
                            0,
                            0 + (i == 0 ? spg::Real(ropeLength * 0.01) / (nParticles - 1) : 0)},
                           i == 0 || i == (nParticles - 1) ? mass / (nParticles - 1) * 0.5 : mass / (nParticles - 1));
    }
    auto springAnchorEnergy = std::make_shared<spg::SpringAnchorEnergy>();

    std::shared_ptr<spg::ParticleGroup::EnergyT> springEnergy;
    if (springType == SpringType::Continuum) {
        springEnergy = std::make_shared<spg::SpringContinuumEnergy>();
    } else if (springType == SpringType::SquaredContinuum) {
        springEnergy = std::make_shared<spg::SpringSquaredContinuumEnergy>();
    } else {
        springEnergy = std::make_shared<spg::SpringEnergy>();
    }
    const spg::Real springAnchorStiffness = 1e3;
    const spg::Real springStiffness = 1e1;
    int anchorParticleIdx = 0;
    springAnchorEnergy->addStencil(
        std::array<int, 1>{anchorParticleIdx}, pGroup.positions0()[anchorParticleIdx], springAnchorStiffness);
    anchorParticleIdx = (nParticles - 1) / 2;
    springAnchorEnergy->addStencil(
        std::array<int, 1>{anchorParticleIdx}, pGroup.positions0()[anchorParticleIdx], springAnchorStiffness);

    auto l_addStencils = [nParticles, ropeLength, springStiffness](auto *energy) {
        for (int i = 1; i < nParticles; i += 1) {
            energy->addStencil(std::array<int, 2>{i - 1, i}, ropeLength / (nParticles - 1), springStiffness);
        }
    };
    if (springType == SpringType::Continuum) {
        l_addStencils(dynamic_cast<spg::SpringContinuumEnergy *>(springEnergy.get()));
    } else if (springType == SpringType::SquaredContinuum) {
        l_addStencils(dynamic_cast<spg::SpringSquaredContinuumEnergy *>(springEnergy.get()));
    } else {
        l_addStencils(dynamic_cast<spg::SpringEnergy *>(springEnergy.get()));
    }
    pGroup.addEnergy(springAnchorEnergy);
    pGroup.addEnergy(springEnergy);
    return pGroup;
}

spg::Real tetrahedralVolume(const spg::Vector3 &x0,
                            const spg::Vector3 &x1,
                            const spg::Vector3 &x2,
                            const spg::Vector3 &x3)
{
    const spg::Vector3 u = x1 - x0;
    const spg::Vector3 v = x2 - x0;
    const spg::Vector3 w = x3 - x0;
    spg::Matrix3 matMatrix;
    matMatrix.col(0) = u;
    matMatrix.col(1) = v;
    matMatrix.col(2) = w;
    return matMatrix.determinant() / 6.0;
}
spg::TetrahedralMesh buildTetrahedralBeamMesh(const float sideX,
                                              const float sideY,
                                              const float sideZ,
                                              const int nVertexX,
                                              const int nVertexY,
                                              const int nVertexZ)
{
    std::vector<spg::Vector3> vertices;
    std::vector<spg::Vector3> vertices0;
    std::vector<spg::Int4> tets;
    vertices.reserve(nVertexX * nVertexY * nVertexZ);
    vertices0.reserve(nVertexX * nVertexY * nVertexZ);
    for (int i = 0; i < nVertexX; ++i) {
        for (int j = 0; j < nVertexY; ++j) {
            for (int k = 0; k < nVertexZ; ++k) {
                vertices.push_back({-sideX * 0.5f + i * sideX / (nVertexX - 1),
                                    -sideY * 0.5f + j * sideY / (nVertexY - 1),
                                    -sideZ * 0.5f + k * sideZ / (nVertexZ - 1)});
                vertices0.push_back({-sideX * 0.5f + i * sideX / (nVertexX - 1),
                                     -sideY * 0.5f + j * sideY / (nVertexY - 1),
                                     -sideZ * 0.5f + k * sideZ / (nVertexZ - 1)});
            }
        }
    }
    // Build hexagonal cells and decompose each in 5 tetrahedra, with mirroring to guarantee continuity
    for (int i = 0; i < nVertexX - 1; ++i) {
        for (int j = 0; j < nVertexY - 1; ++j) {
            for (int k = 0; k < nVertexZ - 1; ++k) {
                int id1 = (i + 0) * nVertexY * nVertexZ + (j + 0) * nVertexZ + k + 0;
                int id2 = (i + 1) * nVertexY * nVertexZ + (j + 0) * nVertexZ + k + 0;
                int id3 = (i + 0) * nVertexY * nVertexZ + (j + 0) * nVertexZ + k + 1;
                int id4 = (i + 1) * nVertexY * nVertexZ + (j + 0) * nVertexZ + k + 1;
                int id5 = (i + 0) * nVertexY * nVertexZ + (j + 1) * nVertexZ + k + 0;
                int id6 = (i + 1) * nVertexY * nVertexZ + (j + 1) * nVertexZ + k + 0;
                int id7 = (i + 0) * nVertexY * nVertexZ + (j + 1) * nVertexZ + k + 1;
                int id8 = (i + 1) * nVertexY * nVertexZ + (j + 1) * nVertexZ + k + 1;

                if (i % 2 == 1) {
                    std::swap(id1, id2);
                    std::swap(id3, id4);
                    std::swap(id5, id6);
                    std::swap(id7, id8);
                }
                if (j % 2 == 1) {
                    std::swap(id1, id5);
                    std::swap(id2, id6);
                    std::swap(id3, id7);
                    std::swap(id4, id8);
                }
                if (k % 2 == 1) {
                    std::swap(id1, id3);
                    std::swap(id2, id4);
                    std::swap(id5, id7);
                    std::swap(id6, id8);
                }

                auto l_positiveVolumeTetIndices = [&vertices0](int id0, int id1, int id2, int id3) -> spg::Int4 {
                    if (tetrahedralVolume(vertices0[id0], vertices0[id1], vertices0[id2], vertices0[id3]) > 0) {
                        return {id0, id1, id2, id3};
                    }
                    return {id1, id0, id2, id3};
                };
                tets.push_back(l_positiveVolumeTetIndices(id1, id3, id4, id7));
                tets.push_back(l_positiveVolumeTetIndices(id7, id8, id4, id6));
                tets.push_back(l_positiveVolumeTetIndices(id4, id2, id1, id6));
                tets.push_back(l_positiveVolumeTetIndices(id1, id5, id7, id6));
                tets.push_back(l_positiveVolumeTetIndices(id7, id4, id6, id1));
            }
        }
    }
    return spg::TetrahedralMesh(vertices, vertices0, tets);
}

spg::ParticleGroup createSpringBeam(const float mass,
                                    const float sideX,
                                    const float sideY,
                                    const float sideZ,
                                    const int nVertexX,
                                    const int nVertexY,
                                    const int nVertexZ,
                                    SpringType springType)
{
    spg::ParticleGroup pGroup;
    spg::TetrahedralMesh tetMesh = buildTetrahedralBeamMesh(sideX, sideY, sideZ, nVertexX, nVertexY, nVertexZ);
    const auto &vertices = tetMesh.vertices();
    const auto &vertices0 = tetMesh.vertices0();
    const auto &edges = tetMesh.edges();
    const auto &tets = tetMesh.tets();
    // Compute lumped particle masses
    std::vector<spg::Real> particleVolumes(vertices.size(), 0);
    spg::Real totalVolume = 0;
    for (const auto &tet : tets) {
        const auto tetVolume =
            tetrahedralVolume(vertices0[tet[0]], vertices0[tet[1]], vertices0[tet[2]], vertices0[tet[3]]);
        totalVolume += tetVolume;
        particleVolumes[tet[0]] += tetVolume / 4.0;
        particleVolumes[tet[1]] += tetVolume / 4.0;
        particleVolumes[tet[2]] += tetVolume / 4.0;
        particleVolumes[tet[3]] += tetVolume / 4.0;
    }
    const spg::Real density = mass / totalVolume;
    for (int i = 0; i < vertices.size(); ++i) {
        pGroup.addParticle(vertices[i], vertices0[i], particleVolumes[i] * density);
    }
    auto springAnchorEnergy = std::make_shared<spg::SpringAnchorEnergy>();
    std::shared_ptr<spg::ParticleGroup::EnergyT> springEnergy;
    if (springType == SpringType::Continuum) {
        springEnergy = std::make_shared<spg::SpringContinuumEnergy>();
    } else if (springType == SpringType::SquaredContinuum) {
        springEnergy = std::make_shared<spg::SpringSquaredContinuumEnergy>();
    } else {
        springEnergy = std::make_shared<spg::SpringEnergy>();
    }

    const spg::Real springStiffness = 1e2;

    // Anchor some nodes
    for (int j = 0; j < nVertexY; ++j) {
        for (int k = 0; k < nVertexZ; ++k) {
            const int id = 0 * nVertexY * nVertexZ + (j + 0) * nVertexZ + k + 0;
            springAnchorEnergy->addStencil(std::array<int, 1>{id}, pGroup.positions0()[id], springStiffness * 1e3);
        }
    }
    auto l_addStencils = [&edges, springStiffness, &pGroup](auto *springEnergy) {
        for (const auto &[id0, id1] : edges) {
            springEnergy->addStencil(std::array<int, 2>{id0, id1},
                                     (pGroup.positions0()[id0] - pGroup.positions0()[id1]).norm(),
                                     springStiffness);
        }
    };
    if (springType == SpringType::Continuum) {
        l_addStencils(dynamic_cast<spg::SpringContinuumEnergy *>(springEnergy.get()));
    } else if (springType == SpringType::SquaredContinuum) {
        l_addStencils(dynamic_cast<spg::SpringSquaredContinuumEnergy *>(springEnergy.get()));
    } else {
        l_addStencils(dynamic_cast<spg::SpringEnergy *>(springEnergy.get()));
    }
    pGroup.addEnergy(springAnchorEnergy);
    pGroup.addEnergy(springEnergy);
    return pGroup;
}

spg::ParticleGroup createTetrahedralBeam(const float mass,
                                         const float sideX,
                                         const float sideY,
                                         const float sideZ,
                                         const int nVertexX,
                                         const int nVertexY,
                                         const int nVertexZ,
                                         const float young,
                                         const float poisson,
                                         FemType femType)
{
    spg::ParticleGroup pGroup;
    spg::TetrahedralMesh tetMesh = buildTetrahedralBeamMesh(sideX, sideY, sideZ, nVertexX, nVertexY, nVertexZ);
    const auto &vertices = tetMesh.vertices();
    const auto &vertices0 = tetMesh.vertices0();
    const auto &tets = tetMesh.tets();
    // Compute lumped particle masses
    std::vector<spg::Real> particleVolumes(vertices.size(), 0);
    spg::Real totalVolume = 0;
    for (const auto &tet : tets) {
        const spg::Vector3 x0{vertices0[tet[0]][0], vertices0[tet[0]][1], vertices0[tet[0]][2]};
        const spg::Vector3 x1{vertices0[tet[1]][0], vertices0[tet[1]][1], vertices0[tet[1]][2]};
        const spg::Vector3 x2{vertices0[tet[2]][0], vertices0[tet[2]][1], vertices0[tet[2]][2]};
        const spg::Vector3 x3{vertices0[tet[3]][0], vertices0[tet[3]][1], vertices0[tet[3]][2]};
        const auto tetVolume =
            tetrahedralVolume(vertices0[tet[0]], vertices0[tet[1]], vertices0[tet[2]], vertices0[tet[3]]);
        totalVolume += tetVolume;
        particleVolumes[tet[0]] += tetVolume / 4.0;
        particleVolumes[tet[1]] += tetVolume / 4.0;
        particleVolumes[tet[2]] += tetVolume / 4.0;
        particleVolumes[tet[3]] += tetVolume / 4.0;
    }
    const spg::Real density = mass / totalVolume;
    for (int i = 0; i < vertices.size(); ++i) {
        pGroup.addParticle(vertices[i], vertices0[i], particleVolumes[i] * density);
    }
    auto springAnchorEnergy = std::make_shared<spg::SpringAnchorEnergy>();
    std::shared_ptr<spg::ParticleGroup::EnergyT> tetEnergy;
    if (femType == FemType::NeoHookean) {
        tetEnergy = std::make_shared<spg::StableNeoHookeanEnergy>();
    } else if (femType == FemType::SquaredNeoHookean) {
        tetEnergy = std::make_shared<spg::StableSquaredNeoHookeanEnergy>();
    } else {
        tetEnergy = std::make_shared<spg::StvkEnergy>();
    }
    const spg::Real springStiffness = 1e4;

    // Anchor some nodes
    for (int j = 0; j < nVertexY; ++j) {
        for (int k = 0; k < nVertexZ; ++k) {
            const int id = 0 * nVertexY * nVertexZ + (j + 0) * nVertexZ + k + 0;
            springAnchorEnergy->addStencil(std::array<int, 1>{id}, pGroup.positions0()[id], springStiffness);
        }
    }
    auto l_addStencils = [&tets, young, poisson, &pGroup](auto *tetEnergy) {
        for (const auto &tetIds : tets) {
            tetEnergy->addStencil(tetIds, young, poisson);
        }
    };
    if (femType == FemType::NeoHookean) {
        l_addStencils(dynamic_cast<spg::StableNeoHookeanEnergy *>(tetEnergy.get()));
    } else if (femType == FemType::SquaredNeoHookean) {
        l_addStencils(dynamic_cast<spg::StableSquaredNeoHookeanEnergy *>(tetEnergy.get()));
    } else {
        l_addStencils(dynamic_cast<spg::StvkEnergy *>(tetEnergy.get()));
    }
    pGroup.addEnergy(springAnchorEnergy);
    pGroup.addEnergy(tetEnergy);
    tetEnergy->preparePrecomputations(pGroup);
    return pGroup;
}

spg::TriangleMesh buildTriangleSquareMesh(const float width,
                                          const float height,
                                          const int nvertexRows,
                                          const int nvertexCols)
{
    std::vector<spg::Vector3> vertices;
    std::vector<spg::Vector3> vertices0;
    std::vector<spg::Int3> faces;
    vertices.reserve(nvertexRows * nvertexCols);
    vertices0.reserve(vertices.size());
    const float rowStep = height / (nvertexRows - 1);
    const float colStep = width / (nvertexCols - 1);
    for (int row = 0; row < nvertexRows; ++row) {
        for (int col = 0; col < nvertexCols; ++col) {
            vertices0.push_back({colStep * col, rowStep * row, 0});
            vertices.push_back({colStep * col, 0, -rowStep * row});
            if (row != 0 && col != 0) {
                if ((row + col) % 2 == 0) {
                    faces.push_back({row * nvertexCols + col - 1,
                                     (row - 1) * nvertexCols + col - 1,
                                     (row - 1) * nvertexCols + col});
                    faces.push_back(
                        {row * nvertexCols + col - 1, (row - 1) * nvertexCols + col, row * nvertexCols + col});
                } else {
                    faces.push_back(
                        {(row - 1) * nvertexCols + col - 1, row * nvertexCols + col, row * nvertexCols + col - 1});
                    faces.push_back(
                        {(row - 1) * nvertexCols + col - 1, (row - 1) * nvertexCols + col, row * nvertexCols + col});
                }
            }
        }
    }
    return spg::TriangleMesh(vertices, vertices0, faces);
}

spg::ParticleGroup createCloth(const float mass,
                               const float width,
                               const float height,
                               const int nvertexRows,
                               const int nvertexCols,
                               MembraneType membraneType,
                               BendingType bendingType)
{
    const auto mesh = buildTriangleSquareMesh(width, height, nvertexRows, nvertexCols);
    spg::ParticleGroup pGroup;
    const auto &vertices = mesh.vertices();
    const auto &vertices0 = mesh.vertices0();
    const auto &faces = mesh.faces();
    const auto &edges = mesh.edges();
    const auto &opposingVerticesPerEdge = mesh.opposingVerticesPerEdge();
    // Compute lumped particle masses
    std::vector<spg::Real> particleAreas(vertices.size(), 0);
    spg::Real totalArea = 0;
    for (const auto &face : faces) {
        const spg::Vector3 x0{vertices0[face[0]][0], vertices0[face[0]][1], vertices0[face[0]][2]};
        const spg::Vector3 x1{vertices0[face[1]][0], vertices0[face[1]][1], vertices0[face[1]][2]};
        const spg::Vector3 x2{vertices0[face[2]][0], vertices0[face[2]][1], vertices0[face[2]][2]};
        const spg::Vector3 e0 = x1 - x0;
        const spg::Vector3 e3 = x2 - x1;
        const auto area = 0.5 * (e0.cross(e3).norm());
        totalArea += area;
        particleAreas[face[0]] += area / 3.0;
        particleAreas[face[1]] += area / 3.0;
        particleAreas[face[2]] += area / 3.0;
    }
    const spg::Real density = mass / totalArea;
    for (int i = 0; i < vertices.size(); ++i) {
        pGroup.addParticle(vertices[i], vertices0[i], particleAreas[i] * density);
    }
    const spg::Real anchorSpringStiffness = 1e5;
    // Anchor corner vertices
    auto springAnchorEnergy = std::make_shared<spg::SpringAnchorEnergy>();
    int anchorParticleIdx = 0;
    springAnchorEnergy->addStencil(
        std::array<int, 1>{anchorParticleIdx}, pGroup.positions()[anchorParticleIdx], anchorSpringStiffness);
    anchorParticleIdx = nvertexRows - 1;
    springAnchorEnergy->addStencil(
        std::array<int, 1>{anchorParticleIdx}, pGroup.positions()[anchorParticleIdx], anchorSpringStiffness);

    // Anchor based on distance to center
    /* auto springAnchorEnergy = std::make_shared<spg::SpringAnchorEnergy>();
    for (int anchorParticleIdx = 0; anchorParticleIdx < vertices.size(); ++anchorParticleIdx) {
        if ((spg::Vector3(
                 vertices[anchorParticleIdx][0], vertices[anchorParticleIdx][1], vertices[anchorParticleIdx][2]) -
             spg::Vector3(width * 0.5, 0, -height * 0.5))
                .norm() < width * 0.25) {
            springAnchorEnergy->addStencil(
                std::array<int, 1>{anchorParticleIdx}, pGroup.positions()[anchorParticleIdx], anchorSpringStiffness);
        }
    } */
    // Anchor first half rows
    /* auto springAnchorEnergy = std::make_shared<spg::SpringAnchorEnergy>();
    for (int anchorParticleIdx = 0; anchorParticleIdx < nvertexCols * static_cast<int>((nvertexRows + 1) * 0.5);
         ++anchorParticleIdx) {
        springAnchorEnergy->addStencil(
            std::array<int, 1>{anchorParticleIdx}, pGroup.positions()[anchorParticleIdx], anchorSpringStiffness);
    } */

    // Add bending energy
    const spg::Real bendingStiffness = 1e-2;
    std::shared_ptr<spg::ParticleGroup::EnergyT> bendingEnergy;
    if (bendingType == BendingType::Discrete) {
        bendingEnergy = std::make_shared<spg::DiscreteBendingEnergy>();
    } else if (bendingType == BendingType::BaraffWitkin) {
        bendingEnergy = std::make_shared<spg::BaraffWitkinBendingEnergy>();
    } else if (bendingType == BendingType::Quadratic) {
        bendingEnergy = std::make_shared<spg::QuadraticBendingEnergy>();
    }

    auto l_addBendingStencil = [&edges, &opposingVerticesPerEdge, bendingStiffness](auto *energy) {
        for (int i = 0; i < edges.size(); ++i) {
            if (opposingVerticesPerEdge[i].size() == 2) {
                energy->addStencil(std::array<int, 4>{edges[i][0],
                                                      edges[i][1],
                                                      opposingVerticesPerEdge[i][0].vertexId,
                                                      opposingVerticesPerEdge[i][1].vertexId},
                                   0.,
                                   bendingStiffness);
            }
        }
    };
    if (bendingType == BendingType::Discrete) {
        l_addBendingStencil(dynamic_cast<spg::DiscreteBendingEnergy *>(bendingEnergy.get()));
    } else if (bendingType == BendingType::BaraffWitkin) {
        l_addBendingStencil(dynamic_cast<spg::BaraffWitkinBendingEnergy *>(bendingEnergy.get()));
    } else if (bendingType == BendingType::Quadratic) {
        l_addBendingStencil(dynamic_cast<spg::QuadraticBendingEnergy *>(bendingEnergy.get()));
    }
    pGroup.addEnergy(bendingEnergy);
    bendingEnergy->preparePrecomputations(pGroup);

    const spg::Real young = 1e2;
    const spg::Real poisson = 0.0;

    std::shared_ptr<spg::ParticleGroup::EnergyT> membraneEnergy;
    if (membraneType == MembraneType::Stvk) {
        membraneEnergy = std::make_shared<spg::MembraneStvkEnergy>();
    } else if (membraneType == MembraneType::BaraffWitkin) {
        membraneEnergy = std::make_shared<spg::MembraneBaraffWitkinEnergy>();
    } else if (membraneType == MembraneType::Choi) {
        membraneEnergy = std::make_shared<spg::MembraneChoiEnergy>();
    }

    auto l_addMembraneStencil = [&faces, young, poisson](auto *energy) {
        for (const auto &face : faces) {
            if constexpr (std::is_same_v<decltype(energy), spg::MembraneBaraffWitkinEnergy *> ||
                          std::is_same_v<decltype(energy), spg::MembraneChoiEnergy *>) {
                energy->addStencil(std::array<int, 3>{face[0], face[1], face[2]}, young, young, young);
            } else {
                energy->addStencil(std::array<int, 3>{face[0], face[1], face[2]}, young, poisson);
            }
        }
    };
    if (membraneType == MembraneType::Stvk) {
        l_addMembraneStencil(dynamic_cast<spg::MembraneStvkEnergy *>(membraneEnergy.get()));
    } else if (membraneType == MembraneType::BaraffWitkin) {
        l_addMembraneStencil(dynamic_cast<spg::MembraneBaraffWitkinEnergy *>(membraneEnergy.get()));
    } else if (membraneType == MembraneType::Choi) {
        l_addMembraneStencil(dynamic_cast<spg::MembraneChoiEnergy *>(membraneEnergy.get()));
    }

    pGroup.addEnergy(membraneEnergy);
    membraneEnergy->preparePrecomputations(pGroup);

    pGroup.addEnergy(springAnchorEnergy);
    return pGroup;
}

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
spg::RigidBodyGroup createAnchoredRigidBody(const spg::Real mass,
                                            const float width,
                                            const float height,
                                            const float depth)
{
    spg::RigidBodyGroup rbGroup;
    auto [vertices, faces] = unitCubeGeometry();
    spg::Matrix3 transform;
    transform.setZero();
    transform(0, 0) = width;
    transform(1, 1) = height;
    transform(2, 2) = depth;
    for (auto &v : vertices) {
        v = transform * v;
    }
    spg::TriangleMesh mesh(vertices, vertices, faces);
    spg::Matrix3 localInertia;
    localInertia.setZero();
    localInertia(0, 0) = 1.0 / 12.0 * mass * (height * height + depth * depth);
    localInertia(1, 1) = 1.0 / 12.0 * mass * (width * width + depth * depth);
    localInertia(2, 2) = 1.0 / 12.0 * mass * (width * width + height * height);

    rbGroup.addBody({0, 0, 0}, {0, 0, 0}, mass, {0, 0, 0}, {0, 0, 0}, localInertia, mesh);

    std::shared_ptr<spg::RBSpringAnchorEnergy> springAnchorEnergy = std::make_shared<spg::RBSpringAnchorEnergy>();
    springAnchorEnergy->addStencil({0}, {0, height * 0.5, 0}, {0, height * 0.5, -1}, 100);
    springAnchorEnergy->addStencil({0}, {width * 0.5, -height * 0.5, 0}, {width, -height * 0.5, 1}, 100);
    springAnchorEnergy->addStencil({0}, {-width * 0.5, -height * 0.5, 0}, {-width, -height * 0.5, 1}, 100);

    rbGroup.addEnergy(springAnchorEnergy);
    return rbGroup;
}

spg::RigidBodyGroup createArticulatedRigidBody(const spg::Real mass,
                                               const float width,
                                               const float height,
                                               const float depth,
                                               const int nBodies)
{
    spg::RigidBodyGroup rbGroup;
    auto [vertices, faces] = unitCubeGeometry();
    spg::Matrix3 transform;
    transform.setZero();
    transform(0, 0) = width;
    transform(1, 1) = height;
    transform(2, 2) = depth;
    for (auto &v : vertices) {
        v = transform * v;
    }
    spg::TriangleMesh mesh(vertices, vertices, faces);
    spg::Matrix3 localInertia;
    localInertia.setZero();
    localInertia(0, 0) = 1.0 / 12.0 * mass * (height * height + depth * depth);
    localInertia(1, 1) = 1.0 / 12.0 * mass * (width * width + depth * depth);
    localInertia(2, 2) = 1.0 / 12.0 * mass * (width * width + height * height);

    rbGroup.addBody({0, 0, 0}, {0, 0, 0}, mass, {0, 0, 0}, {0, 0, 0}, localInertia, mesh);

    std::shared_ptr<spg::RBSpringAnchorEnergy> springAnchorEnergy = std::make_shared<spg::RBSpringAnchorEnergy>();
    springAnchorEnergy->addStencil({0}, {0, height * 0.5, 0}, {0, height * 0.5, 0}, 1e10);
    springAnchorEnergy->addStencil({0}, {0, -height * 0.5, 0}, {0, -height * 0.5, 0}, 1e10);
    springAnchorEnergy->addStencil({0}, {-width * 0.5, 0, 0}, {-width * 0.5, 0, 0}, 1e10);

    std::shared_ptr<spg::RBSpringSquaredEnergy> springEnergy = std::make_shared<spg::RBSpringSquaredEnergy>();
    std::shared_ptr<spg::RBOrientationAlignmentEnergy> alignmentEnergy =
        std::make_shared<spg::RBOrientationAlignmentEnergy>();
    for (int i = 1; i < nBodies; ++i) {
        rbGroup.addBody({width * i, 0, 0}, {width * i, 0, 0}, mass, {0, 0, 0}, {0, 0, 0}, localInertia, mesh);
        springEnergy->addStencil(
            {i - 1, i}, {spg::Vector3{width * 0.5, 0, 0}, spg::Vector3{-width * 0.5, 0, 0}}, 0, 1e10);
        alignmentEnergy->addStencil({i - 1, i}, 1e5);
    }

    rbGroup.addEnergy(springAnchorEnergy);
    if (nBodies > 1) {
        rbGroup.addEnergy(springEnergy);
        rbGroup.addEnergy(alignmentEnergy);
    }
    return rbGroup;
}

spg::RigidBodyGroup createRigidBodyChain(const spg::Real mass,
                                         const float width,
                                         const float height,
                                         const float depth,
                                         const int nBodies)
{
    spg::RigidBodyGroup rbGroup;
    auto [vertices, faces] = unitCubeGeometry();
    spg::Matrix3 transform;
    transform.setZero();
    transform(0, 0) = width;
    transform(1, 1) = height;
    transform(2, 2) = depth;
    for (auto &v : vertices) {
        v = transform * v;
    }
    spg::TriangleMesh mesh(vertices, vertices, faces);
    spg::Matrix3 localInertia;
    localInertia.setZero();
    localInertia(0, 0) = 1.0 / 12.0 * mass * (height * height + depth * depth);
    localInertia(1, 1) = 1.0 / 12.0 * mass * (width * width + depth * depth);
    localInertia(2, 2) = 1.0 / 12.0 * mass * (width * width + height * height);

    rbGroup.addBody({0, 0, 0}, {0, 0, 0}, mass, {0, 0, 0}, {0, 0, 0}, localInertia, mesh);

    std::shared_ptr<spg::RBSpringAnchorEnergy> springAnchorEnergy = std::make_shared<spg::RBSpringAnchorEnergy>();
    springAnchorEnergy->addStencil({0}, {width * 0.5, height * 0.5, 0}, {width * 0.5, height * 0.5, 0}, 100);

    std::shared_ptr<spg::RBSpringEnergy> springEnergy = std::make_shared<spg::RBSpringEnergy>();
    for (int i = 1; i < nBodies; ++i) {
        rbGroup.addBody({-(height + 1) * i, 0, -(height + 1) * i},
                        {-(height + 1) * i, 0, -(height + 1) * i},
                        mass,
                        {0, 0, 0},
                        {0, 0, 0},
                        localInertia,
                        mesh);
        springEnergy->addStencil(
            {i - 1, i},
            {spg::Vector3{-width * 0.5, -height * 0.5, 0}, spg::Vector3{width * 0.5, height * 0.5, 0}},
            std::sqrt(width * 2),
            100);
    }

    rbGroup.addEnergy(springAnchorEnergy);
    if (nBodies > 1) {
        rbGroup.addEnergy(springEnergy);
    }
    return rbGroup;
}

enum class SceneType { Cloth, Rope, SpringBeam, FemBeam, RBChain, AnchoredRB, ArticulatedRB };

std::vector<spg::ParticleGroup> createSceneObjects(SceneType sceneType,
                                                   MembraneType membraneType,
                                                   BendingType bendingType,
                                                   SpringType springType,
                                                   FemType femType,
                                                   int resolutionMultiplier)
{
    float mass = 1.0f;
    if (sceneType == SceneType::Cloth) {
        return {createCloth(mass, 1, 1, 3 * resolutionMultiplier, 3 * resolutionMultiplier, membraneType, bendingType)};
    }
    if (sceneType == SceneType::Rope) {
        return {createRope(mass, 3, 3 * resolutionMultiplier, springType)};
    }
    if (sceneType == SceneType::SpringBeam) {
        return {createSpringBeam(mass,
                                 4,
                                 1,
                                 1,
                                 1 + 4 * resolutionMultiplier,
                                 1 + resolutionMultiplier,
                                 1 + resolutionMultiplier,
                                 springType)};
    }
    if (sceneType == SceneType::FemBeam) {
        return {createTetrahedralBeam(mass,
                                      4,
                                      1,
                                      1,
                                      1 + 4 * resolutionMultiplier,
                                      1 + resolutionMultiplier,
                                      1 + resolutionMultiplier,
                                      1e3F,
                                      0.2F,
                                      femType)};
    }
    throw std::runtime_error("Invalid particle group scene");
}

std::vector<spg::RigidBodyGroup> createRigidBodySceneObjects(SceneType sceneType, int resolutionMultiplier)
{
    float width{1}, height{1}, depth{1}, mass{1};
    if (sceneType == SceneType::RBChain) {
        return {createRigidBodyChain(mass, width, height, depth, 1 * resolutionMultiplier)};
    }
    if (sceneType == SceneType::AnchoredRB) {
        return {createAnchoredRigidBody(mass, width, height, depth)};
    }

    if (sceneType == SceneType::ArticulatedRB) {
        return {createArticulatedRigidBody(mass, width, height, depth, 1 * resolutionMultiplier)};
    }
    throw std::runtime_error("Invalid rigid body group scene");
}

/////////////////////////////////////////
// Polyscope and custom utils for picking
/////////////////////////////////////////

std::shared_ptr<spg::SpringAnchorEnergy> addPickingEnergy(spg::ParticleGroup &pGroup,
                                                          int particleId,
                                                          const spg::Vector3 &pos)
{
    auto springAnchorEnergy = std::make_shared<spg::SpringAnchorEnergy>();
    spg::Real pickingStiffness = 1e4;
    springAnchorEnergy->addStencil(std::array<int, 1>{particleId}, pos, pickingStiffness);
    pGroup.addEnergy(springAnchorEnergy);
    return springAnchorEnergy;
}

std::shared_ptr<spg::RBSpringAnchorEnergy> addPickingEnergy(spg::RigidBodyGroup &rbGroup,
                                                            int bodyId,
                                                            const spg::Vector3 &pos)
{
    auto springAnchorRBEnergy = std::make_shared<spg::RBSpringAnchorEnergy>();
    springAnchorRBEnergy->setName("picking");
    spg::Real pickingStiffness = 100;
    springAnchorRBEnergy->addStencil(
        std::array<int, 1>{bodyId}, rbGroup.globalToLocalRigidBodyFrame(pos, bodyId), pos, pickingStiffness);
    rbGroup.addEnergy(springAnchorRBEnergy);
    return springAnchorRBEnergy;
}

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

int main()
{
    try {
        std::vector<std::shared_ptr<spg::solver::BaseSolver>> solvers;

        // Stuff for picking
        std::unordered_map<const polyscope::Structure *, spg::ParticleGroup *> polyscopeToParticleGroup;
        std::shared_ptr<spg::SpringAnchorEnergy> currentPickingEnergy;
        spg::ParticleGroup *currentPickedParticleGroup = nullptr;
        int currentPickedParticleId = -1;
        std::unordered_map<const polyscope::Structure *, std::pair<spg::RigidBodyGroup *, int>>
            polyscopeToRigidBodyGroupAndBody;
        std::shared_ptr<spg::RBSpringAnchorEnergy> currentRigidBodyPickingEnergy;
        spg::RigidBodyGroup *currentPickedRigidBodyGroup = nullptr;

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
        std::vector<int> solverSubsteps;
        std::vector<int> solverTypes;
        std::vector<const char **> solverTypeNamesLists;
        for (int s = 0; s < solvers.size(); ++s) {
            solverSubsteps.push_back(1);
            solverTypes.push_back(0);
        }

        // Scene selection stuff
        const char *sceneTypeNames[] = {
            "Spring beam", "Cloth", "Rope", "Fem beam", "RB chain", "Anchored RB", "Articulated RB"};
        std::map<int, SceneType> sceneTypeMap;
        sceneTypeMap[0] = SceneType::SpringBeam;
        sceneTypeMap[1] = SceneType::Cloth;
        sceneTypeMap[2] = SceneType::Rope;
        sceneTypeMap[3] = SceneType::FemBeam;
        sceneTypeMap[4] = SceneType::RBChain;
        sceneTypeMap[5] = SceneType::AnchoredRB;
        sceneTypeMap[6] = SceneType::ArticulatedRB;
        int sceneId{0};
        SceneType sceneType = sceneTypeMap[sceneId];

        const char *membraneTypeNames[] = {"StVK", "Baraff Witkin", "Choi"};
        std::map<int, MembraneType> membraneTypeMap;
        membraneTypeMap[0] = MembraneType::Stvk;
        membraneTypeMap[1] = MembraneType::BaraffWitkin;
        membraneTypeMap[2] = MembraneType::Choi;
        const char *bendingTypeNames[] = {"Discrete", "Quadratic", "Baraff Witkin"};
        std::map<int, BendingType> bendingTypeMap;
        bendingTypeMap[0] = BendingType::Discrete;
        bendingTypeMap[1] = BendingType::Quadratic;
        bendingTypeMap[2] = BendingType::BaraffWitkin;
        const char *springTypeNames[] = {"Normal", "Continuum", "Squared Continuum"};
        std::map<int, SpringType> springTypeMap;
        springTypeMap[0] = SpringType::Normal;
        springTypeMap[1] = SpringType::Continuum;
        springTypeMap[2] = SpringType::SquaredContinuum;
        const char *femTypeNames[] = {"StVK", "Neo Hookean", "Squared Neo Hookean"};
        std::map<int, FemType> femTypeMap;
        femTypeMap[0] = FemType::Stvk;
        femTypeMap[1] = FemType::NeoHookean;
        femTypeMap[2] = FemType::SquaredNeoHookean;
        const char *verbosityTypeNames[] = {"None", "Debug", "Performance"};
        std::map<int, spg::Verbosity> verbosityTypeMap;
        verbosityTypeMap[0] = spg::Verbosity::None;
        verbosityTypeMap[1] = spg::Verbosity::Debug;
        verbosityTypeMap[2] = spg::Verbosity::Performance;

        int membraneId{0};
        int bendingId{0};
        int springId{0};
        int femId{0};
        int verbosityId{0};
        MembraneType membraneType = membraneTypeMap[membraneId];
        BendingType bendingType = bendingTypeMap[bendingId];
        SpringType springType = springTypeMap[springId];
        FemType femType = femTypeMap[femId];
        spg::Verbosity verbosity{verbosityTypeMap[verbosityId]};

        int simObjectResolutionMultiplier{1};

        auto l_callback = [&]() {
            auto l_updateView = [&solvers]() {
                for (int solverId = 0; solverId < solvers.size(); ++solverId) {
                    const auto &solver = solvers[solverId];
                    for (int i = 0; i < solver->particleGroups().size(); ++i) {
                        auto *pc = polyscope::getPointCloud("solver" + std::to_string(solverId) + "simObj" +
                                                            std::to_string(i));
                        pc->updatePointPositions(solver->particleGroups()[i].positions());
                        if (polyscope::hasCurveNetwork("solver" + std::to_string(solverId) + "simObj" +
                                                       std::to_string(i) + "springs")) {
                            auto *cn = polyscope::getCurveNetwork("solver" + std::to_string(solverId) + "simObj" +
                                                                  std::to_string(i) + "springs");
                            cn->updateNodePositions(solver->particleGroups()[i].positions());
                        }
                        if (polyscope::hasSurfaceMesh("solver" + std::to_string(solverId) + "simObj" +
                                                      std::to_string(i) + "membrane")) {
                            auto *sm = polyscope::getSurfaceMesh("solver" + std::to_string(solverId) + "simObj" +
                                                                 std::to_string(i) + "membrane");
                            sm->updateVertexPositions(solver->particleGroups()[i].positions());
                        }
                        if (polyscope::hasVolumeMesh("solver" + std::to_string(solverId) + "simObj" +
                                                     std::to_string(i) + "fem")) {
                            auto *vm = polyscope::getVolumeMesh("solver" + std::to_string(solverId) + "simObj" +
                                                                std::to_string(i) + "fem");
                            vm->updateVertexPositions(solver->particleGroups()[i].positions());
                        }
                    }
                    for (int i = 0; i < solver->rbGroups().size(); ++i) {
                        auto &rbGroup = solver->rbGroups()[i];
                        for (int j = 0; j < rbGroup.size(); ++j) {
                            auto *sm = polyscope::getSurfaceMesh("solver" + std::to_string(solverId) + "rbGroup" +
                                                                 std::to_string(i) + "rb" + std::to_string(j));
                            const auto &pos = rbGroup.positions()[j];
                            const auto &rotation = rbGroup.rotationMatrices()[j];
                            glm::mat4 mat(rotation(0, 0),
                                          rotation(0, 1),
                                          rotation(0, 2),
                                          pos(0),
                                          rotation(1, 0),
                                          rotation(1, 1),
                                          rotation(1, 2),
                                          pos(1),
                                          rotation(2, 0),
                                          rotation(2, 1),
                                          rotation(2, 2),
                                          pos(2),
                                          0,
                                          0,
                                          0,
                                          1);
                            mat = glm::transpose(mat);
                            sm->setTransform(mat);
                        }
                        for (auto energy : rbGroup.energies()) {
                            if (auto springEnergy = dynamic_cast<spg::RBSpringEnergy *>(energy.get());
                                springEnergy != nullptr) {
                                std::vector<spg::Vector3> springPositions;
                                for (int s = 0; s < springEnergy->nStencils(); ++s) {
                                    auto rbIds = springEnergy->stencils()[s];
                                    springPositions.push_back(rbGroup.positions()[rbIds[0]] +
                                                              rbGroup.rotationMatrices()[rbIds[0]] *
                                                                  (springEnergy->localRBPoints()[s][0]));
                                    springPositions.push_back(rbGroup.positions()[rbIds[1]] +
                                                              rbGroup.rotationMatrices()[rbIds[1]] *
                                                                  (springEnergy->localRBPoints()[s][1]));
                                }
                                auto *cn = polyscope::getCurveNetwork("solver" + std::to_string(solverId) + "rbGroup" +
                                                                      std::to_string(i) + "springs");
                                cn->updateNodePositions(springPositions);
                            } else if (auto springEnergy = dynamic_cast<spg::RBSpringSquaredEnergy *>(energy.get());
                                       springEnergy != nullptr) {
                                std::vector<spg::Vector3> springPositions;
                                for (int s = 0; s < springEnergy->nStencils(); ++s) {
                                    auto rbIds = springEnergy->stencils()[s];
                                    springPositions.push_back(rbGroup.positions()[rbIds[0]] +
                                                              rbGroup.rotationMatrices()[rbIds[0]] *
                                                                  (springEnergy->localRBPoints()[s][0]));
                                    springPositions.push_back(rbGroup.positions()[rbIds[1]] +
                                                              rbGroup.rotationMatrices()[rbIds[1]] *
                                                                  (springEnergy->localRBPoints()[s][1]));
                                }
                                auto *cn = polyscope::getCurveNetwork("solver" + std::to_string(solverId) + "rbGroup" +
                                                                      std::to_string(i) + "squaredSprings");
                                cn->updateNodePositions(springPositions);
                            } else if (auto springAnchorEnergy =
                                           dynamic_cast<spg::RBSpringAnchorEnergy *>(energy.get());
                                       springAnchorEnergy != nullptr) {
                                std::vector<spg::Vector3> springPositions;
                                for (int s = 0; s < springAnchorEnergy->nStencils(); ++s) {
                                    auto rbIds = springAnchorEnergy->stencils()[s];
                                    springPositions.push_back(rbGroup.positions()[rbIds[0]] +
                                                              rbGroup.rotationMatrices()[rbIds[0]] *
                                                                  (springAnchorEnergy->localRBPoints()[s]));
                                    springPositions.push_back(springAnchorEnergy->anchors()[s]);
                                }
                                auto *cn = polyscope::getCurveNetwork("solver" + std::to_string(solverId) + "rbGroup" +
                                                                      std::to_string(i) + "anchor-springs-" +
                                                                      springAnchorEnergy->name());
                                cn->updateNodePositions(springPositions);
                            }
                        }
                    }
                }
            };
            auto l_registerSolver = [&polyscopeToParticleGroup, &polyscopeToRigidBodyGroupAndBody](auto &solver,
                                                                                                   int solverId) {
                for (int i = 0; i < solver->particleGroups().size(); ++i) {
                    auto &pGroup = solver->particleGroups()[i];
                    polyscope::PointCloud *pc = polyscope::registerPointCloud(
                        "solver" + std::to_string(solverId) + "simObj" + std::to_string(i), pGroup.positions());
                    pc->setPointRadius(0.01);
                    polyscopeToParticleGroup[pc] = &pGroup;

                    /* std::vector<spg::coloring::FlatStencils> flatStencilsSet;
                    for (const auto &energy : pGroup.energies()) {
                        flatStencilsSet.push_back({energy->flatStencils(), energy->stencilSize()});
                    }
                    auto vertexColors = spg::coloring::colorVertices(pGroup.nParticles(), flatStencilsSet);
                    pc->addScalarQuantity("color", vertexColors); */

                    auto l_registerCurveNetwork = [&pGroup, solverId, i](const auto *energy) {
                        polyscope::CurveNetwork *cn = polyscope::registerCurveNetwork(
                            "solver" + std::to_string(solverId) + "simObj" + std::to_string(i) + "springs",
                            pGroup.positions(),
                            energy->stencils());
                        /* auto colors = spg::coloring::colorStencils(
                            {springEnergy->flatStencils(), springEnergy->stencilSize()});
                        cn->addEdgeScalarQuantity("color", colors); */
                    };
                    auto l_registerSurfaceMesh = [&pGroup, solverId, i](const auto *energy) {
                        polyscope::SurfaceMesh *sm = polyscope::registerSurfaceMesh(
                            "solver" + std::to_string(solverId) + "simObj" + std::to_string(i) + "membrane",
                            pGroup.positions(),
                            energy->stencils());
                        /* auto colors = spg::coloring::colorStencils(
                            {membraneEnergy->flatStencils(), membraneEnergy->stencilSize()});
                        sm->addFaceScalarQuantity("color", colors); */
                    };
                    auto l_registerVolumeMesh = [&pGroup, solverId, i](const auto *energy) {
                        polyscope::VolumeMesh *vm = polyscope::registerTetMesh(
                            "solver" + std::to_string(solverId) + "simObj" + std::to_string(i) + "fem",
                            pGroup.positions(),
                            energy->stencils());
                        /* auto colors =
                            spg::coloring::colorStencils({femEnergy->flatStencils(), femEnergy->stencilSize()});
                        vm->addCellScalarQuantity("color", colors); */
                    };
                    for (auto energy : pGroup.energies()) {
                        if (auto springEnergy = dynamic_cast<spg::SpringEnergy *>(energy.get());
                            springEnergy != nullptr) {
                            l_registerCurveNetwork(springEnergy);
                        } else if (auto springContinuumEnergy =
                                       dynamic_cast<spg::SpringContinuumEnergy *>(energy.get());
                                   springContinuumEnergy != nullptr) {
                            l_registerCurveNetwork(springContinuumEnergy);
                        } else if (auto springSquaredContinuumEnergy =
                                       dynamic_cast<spg::SpringSquaredContinuumEnergy *>(energy.get());
                                   springSquaredContinuumEnergy != nullptr) {
                            l_registerCurveNetwork(springSquaredContinuumEnergy);
                        } else if (auto membraneStvkEnergy = dynamic_cast<spg::MembraneStvkEnergy *>(energy.get());
                                   membraneStvkEnergy != nullptr) {
                            l_registerSurfaceMesh(membraneStvkEnergy);
                        } else if (auto membraneChoiEnergy = dynamic_cast<spg::MembraneChoiEnergy *>(energy.get());
                                   membraneChoiEnergy != nullptr) {
                            l_registerSurfaceMesh(membraneChoiEnergy);
                        } else if (auto membraneBaraffWitkinEnergy =
                                       dynamic_cast<spg::MembraneBaraffWitkinEnergy *>(energy.get());
                                   membraneBaraffWitkinEnergy != nullptr) {
                            l_registerSurfaceMesh(membraneBaraffWitkinEnergy);
                        } else if (auto stableNeoHookeanEnergy =
                                       dynamic_cast<spg::StableNeoHookeanEnergy *>(energy.get());
                                   stableNeoHookeanEnergy != nullptr) {
                            l_registerVolumeMesh(stableNeoHookeanEnergy);
                        } else if (auto StableSquaredNeoHookeanEnergy =
                                       dynamic_cast<spg::StableSquaredNeoHookeanEnergy *>(energy.get());
                                   StableSquaredNeoHookeanEnergy != nullptr) {
                            l_registerVolumeMesh(StableSquaredNeoHookeanEnergy);
                        } else if (auto stvkEnergy = dynamic_cast<spg::StvkEnergy *>(energy.get());
                                   stvkEnergy != nullptr) {
                            l_registerVolumeMesh(stvkEnergy);
                        }
                    }
                }
                for (int i = 0; i < solver->rbGroups().size(); ++i) {
                    spg::RigidBodyGroup &rbGroup = solver->rbGroups()[i];
                    for (int j = 0; j < rbGroup.size(); ++j) {
                        polyscope::SurfaceMesh *sm =
                            polyscope::registerSurfaceMesh("solver" + std::to_string(solverId) + "rbGroup" +
                                                               std::to_string(i) + "rb" + std::to_string(j),
                                                           rbGroup.visualMeshes()[j].vertices(),
                                                           rbGroup.visualMeshes()[j].faces());
                        polyscopeToRigidBodyGroupAndBody[sm] = std::make_pair(&rbGroup, j);
                        const auto &pos = rbGroup.positions()[i];
                        const auto &rotation = rbGroup.rotationMatrices()[i];
                        glm::mat4 mat(rotation(0, 0),
                                      rotation(0, 1),
                                      rotation(0, 2),
                                      pos(0),
                                      rotation(1, 0),
                                      rotation(1, 1),
                                      rotation(1, 2),
                                      pos(1),
                                      rotation(2, 0),
                                      rotation(2, 1),
                                      rotation(2, 2),
                                      pos(2),
                                      0,
                                      0,
                                      0,
                                      1);
                        mat = glm::transpose(mat);
                        sm->setTransform(mat);
                    }
                    for (auto energy : rbGroup.energies()) {
                        if (auto springEnergy = dynamic_cast<spg::RBSpringEnergy *>(energy.get());
                            springEnergy != nullptr) {
                            std::vector<spg::Vector3> springPositions;
                            std::vector<std::array<int, 2>> springIndices;
                            for (int s = 0; s < springEnergy->nStencils(); ++s) {
                                auto rbIds = springEnergy->stencils()[s];
                                springPositions.push_back(rbGroup.positions()[rbIds[0]] +
                                                          rbGroup.rotationMatrices()[rbIds[0]] *
                                                              (springEnergy->localRBPoints()[s][0]));
                                springPositions.push_back(rbGroup.positions()[rbIds[1]] +
                                                          rbGroup.rotationMatrices()[rbIds[1]] *
                                                              (springEnergy->localRBPoints()[s][1]));
                                springIndices.push_back({s * 2, s * 2 + 1});
                            }
                            polyscope::registerCurveNetwork(
                                "solver" + std::to_string(solverId) + "rbGroup" + std::to_string(i) + "springs",
                                springPositions,
                                springIndices);
                        } else if (auto springEnergy = dynamic_cast<spg::RBSpringSquaredEnergy *>(energy.get());
                                   springEnergy != nullptr) {
                            std::vector<spg::Vector3> springPositions;
                            std::vector<std::array<int, 2>> springIndices;
                            for (int s = 0; s < springEnergy->nStencils(); ++s) {
                                auto rbIds = springEnergy->stencils()[s];
                                springPositions.push_back(rbGroup.positions()[rbIds[0]] +
                                                          rbGroup.rotationMatrices()[rbIds[0]] *
                                                              (springEnergy->localRBPoints()[s][0]));
                                springPositions.push_back(rbGroup.positions()[rbIds[1]] +
                                                          rbGroup.rotationMatrices()[rbIds[1]] *
                                                              (springEnergy->localRBPoints()[s][1]));
                                springIndices.push_back({s * 2, s * 2 + 1});
                            }
                            polyscope::registerCurveNetwork(
                                "solver" + std::to_string(solverId) + "rbGroup" + std::to_string(i) + "squaredSprings",
                                springPositions,
                                springIndices);
                        } else if (auto springAnchorEnergy = dynamic_cast<spg::RBSpringAnchorEnergy *>(energy.get());
                                   springAnchorEnergy != nullptr) {
                            std::vector<spg::Vector3> springPositions;
                            std::vector<std::array<int, 2>> springIndices;
                            for (int s = 0; s < springAnchorEnergy->nStencils(); ++s) {
                                auto rbIds = springAnchorEnergy->stencils()[s];
                                springPositions.push_back(rbGroup.positions()[rbIds[0]] +
                                                          rbGroup.rotationMatrices()[rbIds[0]] *
                                                              (springAnchorEnergy->localRBPoints()[s]));
                                springPositions.push_back(springAnchorEnergy->anchors()[s]);
                                springIndices.push_back({s * 2, s * 2 + 1});
                            }
                            polyscope::registerCurveNetwork("solver" + std::to_string(solverId) + "rbGroup" +
                                                                std::to_string(i) + "anchor-springs-" +
                                                                springAnchorEnergy->name(),
                                                            springPositions,
                                                            springIndices);
                        }
                    }
                }
            };
            auto l_unregisterSolver = [](const auto &solver, int solverId) {
                for (int i = 0; i < solver->particleGroups().size(); ++i) {
                    polyscope::removePointCloud("solver" + std::to_string(solverId) + "simObj" + std::to_string(i));
                    polyscope::removeCurveNetwork("solver" + std::to_string(solverId) + "simObj" + std::to_string(i) +
                                                  "springs");
                    polyscope::removeSurfaceMesh("solver" + std::to_string(solverId) + "simObj" + std::to_string(i) +
                                                 "membrane");
                    polyscope::removeVolumeMesh("solver" + std::to_string(solverId) + "simObj" + std::to_string(i) +
                                                "fem");
                }
                for (int i = 0; i < solver->rbGroups().size(); ++i) {
                    auto &rbGroup = solver->rbGroups()[i];
                    for (int j = 0; j < rbGroup.size(); ++j) {
                        polyscope::removeSurfaceMesh("solver" + std::to_string(solverId) + "rbGroup" +
                                                     std::to_string(i) + "rb" + std::to_string(j));
                        polyscope::removeCurveNetwork("solver" + std::to_string(solverId) + "rbGroup" +
                                                      std::to_string(j) + "springs");
                        for (auto energy : rbGroup.energies()) {
                            if (auto springAnchorEnergy = dynamic_cast<spg::RBSpringAnchorEnergy *>(energy.get());
                                springAnchorEnergy != nullptr) {
                                polyscope::removeCurveNetwork("solver" + std::to_string(solverId) + "rbGroup" +
                                                              std::to_string(j) + "anchor-springs-" +
                                                              springAnchorEnergy->name());
                            }
                        }
                    }
                }
            };
            auto l_addSolver = [&]() {
                auto &solver = solvers.emplace_back(std::make_shared<spg::solver::XPBD>());
                solverTypes.push_back(0);
                solver->setNumSubsteps(1);
                solver->setVerbosity(verbosity);
                solverSubsteps.push_back(1);
                if (sceneType != SceneType::RBChain && sceneType != SceneType::AnchoredRB) {
                    const auto objects = createSceneObjects(
                        sceneType, membraneType, bendingType, springType, femType, simObjectResolutionMultiplier);
                    for (const auto &object : objects) {
                        solver->addObject(object);
                    }
                } else {
                    const auto objects = createRigidBodySceneObjects(sceneType, simObjectResolutionMultiplier);
                    for (const auto &object : objects) {
                        solver->addObject(object);
                    }
                }
                solver->setDt(dt);
                l_registerSolver(solver, static_cast<int>(solvers.size()) - 1);
            };
            auto l_resetSolver = [&](const int solverId) {
                const auto solverType = solverTypes[solverId];
                auto &solver = solvers[solverId];
                int typeIndex = 0;
                if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::XPBD>();
                } else if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::XPBD>(false);
                } else if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::VBD>();
                } else if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::ImplicitEulerBaraffWitkin>();
                } else if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::BDF2>();
                } else if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::ImplicitEulerNewtonRobust>();
                } else if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::ImplicitEulerNewtonDx>();
                } else if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::ImplicitEulerNewtonDv>();
                } else if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::SimplecticEuler>();
                } else if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::StaticNewton>();
                } else if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::QuasiStaticNewtonRobust>();
                } else if (solverType == typeIndex++) {
                    solver = std::make_shared<spg::solver::QuasiStaticNewton>();
                }
                solver->setNumSubsteps(solverSubsteps[solverId]);
                solver->setVerbosity(verbosity);
                if (sceneType != SceneType::RBChain && sceneType != SceneType::AnchoredRB &&
                    sceneType != SceneType::ArticulatedRB) {
                    const auto objects = createSceneObjects(
                        sceneType, membraneType, bendingType, springType, femType, simObjectResolutionMultiplier);
                    for (const auto &object : objects) {
                        solver->addObject(object);
                    }
                } else {
                    const auto objects = createRigidBodySceneObjects(sceneType, simObjectResolutionMultiplier);
                    for (const auto &object : objects) {
                        solver->addObject(object);
                    }
                }
                solver->setDt(dt);
            };
            auto l_resetAllSolvers = [&]() {
                for (int i = 0; i < solvers.size(); ++i) {
                    l_unregisterSolver(solvers[i], i);
                    l_resetSolver(i);
                    l_registerSolver(solvers[i], i);
                }
                l_updateView();
            };
            auto l_removeSolver = [&solvers, &solverTypes, &solverSubsteps, dt, l_unregisterSolver]() {
                l_unregisterSolver(solvers.back(), static_cast<int>(solvers.size()) - 1);
                solvers.pop_back();
                solverTypes.pop_back();
                solverSubsteps.pop_back();
            };
            if (solvers.empty()) {
                l_addSolver();
            }
            ImGui::PushItemWidth(200);
            if (ImGui::Button("Add solver")) {
                l_addSolver();
            }
            if (ImGui::Button("Remove solver")) {
                if (solvers.size() > 1) {
                    l_removeSolver();
                }
            }
            ImGui::Checkbox("Simulate [Space]", &run);

            if (ImGui::Combo(std::string("Verbosity").c_str(),
                             &verbosityId,
                             verbosityTypeNames,
                             IM_ARRAYSIZE(verbosityTypeNames))) {
                verbosity = verbosityTypeMap[verbosityId];
                for (auto &solver : solvers) {
                    solver->setVerbosity(verbosity);
                }
            }
            if (ImGui::Combo(std::string("Scene").c_str(), &sceneId, sceneTypeNames, IM_ARRAYSIZE(sceneTypeNames))) {
                sceneType = sceneTypeMap[sceneId];
                simObjectResolutionMultiplier = 1;
                l_resetAllSolvers();
                polyscope::view::resetCameraToHomeView();
            }
            if (ImGui::InputInt("Obj resolution multiplier", &simObjectResolutionMultiplier)) {
                simObjectResolutionMultiplier = std::max(1, simObjectResolutionMultiplier);
                l_resetAllSolvers();
            }
            if (ImGui::Combo(std::string("Membrane model").c_str(),
                             &membraneId,
                             membraneTypeNames,
                             IM_ARRAYSIZE(membraneTypeNames))) {
                membraneType = membraneTypeMap[membraneId];
                l_resetAllSolvers();
            }
            if (ImGui::Combo(std::string("Bending model").c_str(),
                             &bendingId,
                             bendingTypeNames,
                             IM_ARRAYSIZE(bendingTypeNames))) {
                bendingType = bendingTypeMap[bendingId];
                l_resetAllSolvers();
            }
            if (ImGui::Combo(
                    std::string("Spring model").c_str(), &springId, springTypeNames, IM_ARRAYSIZE(springTypeNames))) {
                springType = springTypeMap[springId];
                l_resetAllSolvers();
            }
            if (ImGui::Combo(std::string("Fem model").c_str(), &femId, femTypeNames, IM_ARRAYSIZE(femTypeNames))) {
                femType = femTypeMap[femId];
                l_resetAllSolvers();
            }
            ImGui::Dummy(ImVec2(0.0f, 20.0f));
            const char *solverTypeNames[] = {"XPBD parallel Gauss-Seidel",
                                             "XPBD serial Gauss-Seidel",
                                             "Vertex Block Descent",
                                             "Baraff-Witkin",
                                             "BDF2",
                                             "Implicit Euler Newton robust",
                                             "Implicit Euler Newton on x",
                                             "Implicit Euler Newton on v",
                                             "Simplectic Euler",
                                             "Static Newton",
                                             "Quasi-static Newton robust",
                                             "Quasi-static Newton"};
            for (int s = 0; s < solvers.size(); ++s) {
                auto &solverType = solverTypes[s];
                if (ImGui::Combo(std::string("Solver " + std::to_string(s) + " type").c_str(),
                                 &solverType,
                                 solverTypeNames,
                                 IM_ARRAYSIZE(solverTypeNames))) {
                    l_unregisterSolver(solvers[s], s);
                    l_resetSolver(s);
                    l_registerSolver(solvers[s], s);
                }
                l_updateView();
            }
            for (int i = 0; i < solverSubsteps.size(); ++i) {
                auto &substeps = solverSubsteps[i];
                if (ImGui::InputInt(std::string("Solver " + std::to_string(i) + " substeps").c_str(), &substeps)) {
                    solvers[i]->setNumSubsteps(substeps);
                }
            }
            if (ImGui::InputFloat("dt", &dt, 0, 0, "%.6f")) {
                dt = std::max(dt, 1e-6f);
                for (auto &solver : solvers) {
                    solver->setDt(dt);
                }
            }
            if (ImGui::Button("Double mass")) {
                for (auto &solver : solvers) {
                    for (auto &obj : solver->particleGroups()) {
                        obj.scaleMasses(2);
                    }
                    for (auto &obj : solver->rbGroups()) {
                        obj.scaleMasses(2);
                    }
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Halve mass")) {
                for (auto &solver : solvers) {
                    for (auto &obj : solver->particleGroups()) {
                        obj.scaleMasses(0.5);
                    }
                    for (auto &obj : solver->rbGroups()) {
                        obj.scaleMasses(0.5);
                    }
                }
            }
            if (ImGui::Button("Reset solvers [r]")) {
                for (auto &solver : solvers) {
                    solver->reset();
                }
                l_updateView();
            }

            if (ImGui::Button("Single Step [s]")) {
                if (!run) {
                    for (auto &solver : solvers) {
                        solver->step();
                    }
                    l_updateView();
                }
            }
            if (ImGui::Button("Copy positions from first solver to the rest")) {
                // Assume there is a single simObject per solver, just to simplify
                // Note: If the solvers require additional info to work properly (e.g., BDF2 requires info about the
                // previous step), just copying the current state can lead to weird behaviors
                if (!solvers.front()->particleGroups().empty()) {
                    const auto &positions = solvers.front()->particleGroups().front().positions();
                    const auto &velocities = solvers.front()->particleGroups().front().velocities();
                    for (auto &solver : solvers) {
                        solver->particleGroups().front().positions() = positions;
                        solver->particleGroups().front().velocities() = velocities;
                    }
                }
                if (!solvers.front()->rbGroups().empty()) {
                    const auto &positions = solvers.front()->rbGroups().front().positions();
                    const auto &velocities = solvers.front()->rbGroups().front().velocities();
                    const auto &thetas = solvers.front()->rbGroups().front().thetas();
                    const auto &omegas = solvers.front()->rbGroups().front().omegas();
                    for (auto &solver : solvers) {
                        solver->rbGroups().front().positions() = positions;
                        solver->rbGroups().front().velocities() = velocities;
                        solver->rbGroups().front().thetas() = thetas;
                        solver->rbGroups().front().omegas() = omegas;
                        solver->rbGroups().front().updateRotationMatrices();
                        solver->rbGroups().front().updateInertias();
                    }
                }
                if (!run) {
                    l_updateView();
                }
            }
            if (run) {
                for (auto &solver : solvers) {
                    solver->step();
                }
                l_updateView();
            }

            if (ImGui::Button("Center camera")) {
                polyscope::view::resetCameraToHomeView();
            }

            ImGui::Dummy(ImVec2(0.0f, 20.0f));
            ImGui::Text("Mid-click and drag on a particle or a rigid body to pull it");
            ImGui::Dummy(ImVec2(0.0f, 20.0f));

            ImGui::PopItemWidth();

            if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_Space)) {
                run = !run;
            }

            if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_S)) {
                if (!run) {
                    for (auto &solver : solvers) {
                        solver->step();
                    }
                    l_updateView();
                }
            }

            if (ImGui::IsKeyPressed(ImGuiKey::ImGuiKey_R)) {
                for (auto &solver : solvers) {
                    solver->reset();
                }
                l_updateView();
            }

            // Dirty code to flag specific solver types if a topology change occurs to recompute required internal
            // structures
            auto l_flagSolversForTopologyChanges = [&]() {
                for (int i = 0; i < solvers.size(); ++i) {
                    if (auto *solver = dynamic_cast<spg::solver::VBD *>(solvers[i].get()); solver != nullptr) {
                        solver->requirePrecomputationUpdate();
                    }
                }
            };
            ImGuiIO &io = ImGui::GetIO();
            if (io.MouseClicked[2]) {  // if the center mouse button was clicked
                // gather values
                glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
                std::pair<polyscope::Structure *, size_t> pickPair = polyscope::pick::evaluatePickQuery(
                    static_cast<int>(screenCoords.x), static_cast<int>(screenCoords.y));
                if (pickPair.first != nullptr &&
                    polyscopeToParticleGroup.find(pickPair.first) != polyscopeToParticleGroup.end()) {
                    const auto pos = screenCoordsToWorldPosition(screenCoords, currentPickedDepth);
                    // Not sure why it can return 1 when picking something, but it is the case, ignore the hit
                    if (currentPickedDepth != 1.) {
                        currentPickedParticleGroup = polyscopeToParticleGroup[pickPair.first];
                        currentPickedParticleId = static_cast<int>(pickPair.second);
                        currentPickedPos = {pos.x, pos.y, pos.z};
                        pickAnchorPos = currentPickedParticleGroup->positions()[pickPair.second];
                        currentPickingEnergy = addPickingEnergy(
                            *currentPickedParticleGroup, static_cast<int>(pickPair.second), pickAnchorPos);
                        l_flagSolversForTopologyChanges();
                    }
                } else if (pickPair.first != nullptr && polyscopeToRigidBodyGroupAndBody.find(pickPair.first) !=
                                                            polyscopeToRigidBodyGroupAndBody.end()) {
                    const auto pos = screenCoordsToWorldPosition(screenCoords, currentPickedDepth);
                    // Not sure why it can return 1 when picking something, but it is the case, ignore the hit
                    if (currentPickedDepth != 1.) {
                        auto [rbGroup, bodyId] = polyscopeToRigidBodyGroupAndBody[pickPair.first];
                        currentPickedPos = {pos.x, pos.y, pos.z};
                        pickAnchorPos = currentPickedPos;
                        currentPickedRigidBodyGroup = rbGroup;
                        for (int i = 0; i < solvers.size(); ++i) {
                            l_unregisterSolver(solvers[i], i);
                        }
                        currentRigidBodyPickingEnergy =
                            addPickingEnergy(*currentPickedRigidBodyGroup, bodyId, pickAnchorPos);
                        for (int i = 0; i < solvers.size(); ++i) {
                            l_registerSolver(solvers[i], i);
                        }
                        l_updateView();
                        l_flagSolversForTopologyChanges();
                    }
                }
            } else if (io.MouseReleased[2]) {
                currentPickedDepth = -1;
                if (currentPickedParticleGroup != nullptr) {
                    currentPickedParticleGroup->removeEnergy(currentPickingEnergy);
                    currentPickedParticleGroup = nullptr;
                    currentPickingEnergy = nullptr;
                    l_flagSolversForTopologyChanges();
                } else if (currentPickedRigidBodyGroup != nullptr) {
                    for (int i = 0; i < solvers.size(); ++i) {
                        l_unregisterSolver(solvers[i], i);
                    }
                    currentPickedRigidBodyGroup->removeEnergy(currentRigidBodyPickingEnergy);
                    for (int i = 0; i < solvers.size(); ++i) {
                        l_registerSolver(solvers[i], i);
                    }
                    l_updateView();
                    currentPickedRigidBodyGroup = nullptr;
                    currentRigidBodyPickingEnergy = nullptr;
                    l_flagSolversForTopologyChanges();
                }
            } else if (io.MouseDown[2] && (io.MouseDelta.x != 0 || io.MouseDelta.y != 0) && currentPickedDepth > 0) {
                glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
                glm::vec3 worldPos = screenCoordsToWorldPosition(screenCoords, currentPickedDepth);

                if (currentPickedParticleGroup != nullptr) {
                    spg::Vector3 disp = spg::Vector3(worldPos.x, worldPos.y, worldPos.z) - currentPickedPos;
                    currentPickedPos = {worldPos.x, worldPos.y, worldPos.z};
                    pickAnchorPos += disp;
                    currentPickingEnergy->updateAnchor(0, pickAnchorPos);
                } else if (currentPickedRigidBodyGroup != nullptr) {
                    spg::Vector3 disp = spg::Vector3(worldPos.x, worldPos.y, worldPos.z) - currentPickedPos;
                    currentPickedPos = {worldPos.x, worldPos.y, worldPos.z};
                    pickAnchorPos += disp;
                    currentRigidBodyPickingEnergy->updateAnchor(0, pickAnchorPos);
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