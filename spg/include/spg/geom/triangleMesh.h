#pragma once

#include <spg/types.h>
#include <stdexcept>

namespace spg
{
class TriangleMesh
{
public:
    struct OpposingVertex {
        int triangleId;
        int vertexId;
    };

    TriangleMesh(const std::vector<Vector3> &vertices,
                 const std::vector<Vector3> &vertices0,
                 const std::vector<Int3> &faces)
        : m_vertices(vertices)
        , m_vertices0(vertices0)
        , m_faces(faces)
    {
        const int nVertices{static_cast<int>(m_vertices.size())};
        if (std::find_if(m_faces.begin(), m_faces.end(), [nVertices](const auto &face) {
                return std::any_of(face.begin(), face.end(), [nVertices](const auto &index) {
                    return index < 0 || index >= nVertices;
                });
            }) != m_faces.end()) {
            throw std::runtime_error("Triangle mesh face contains index to invalid vertex");
        }
        if (std::find_if(m_faces.begin(), m_faces.end(), [nVertices](const auto &face) {
                return face[0] == face[1] || face[1] == face[2] || face[2] == face[0];
            }) != m_faces.end()) {
            throw std::runtime_error("Triangle mesh face contains repeated vertex index");
        }
        // Compute edges info
        std::map<Int2, std::vector<OpposingVertex>> edges;
        for (int f = 0; f < static_cast<int>(m_faces.size()); ++f) {
            const auto &face{m_faces[f]};
            for (int i = 0; i < 3; ++i) {
                const int id0 = face[i];
                const int id1 = face[(i + 1) % 3];
                const int id2 = face[(i + 2) % 3];
                edges[{std::min(id0, id1), std::max(id0, id1)}].push_back({f, id2});
            }
        }
        for (const auto &[edge, opposingVertices] : edges) {
            m_edges.emplace_back(edge);
            m_opposingVerticesPerEdge.emplace_back(opposingVertices);
        }
    }
    const std::vector<Vector3> &vertices() const { return m_vertices; }
    const std::vector<Vector3> &vertices0() const { return m_vertices0; }
    const std::vector<Int3> &faces() const { return m_faces; }
    const std::vector<Int2> &edges() const { return m_edges; }
    const std::vector<std::vector<OpposingVertex>> &opposingVerticesPerEdge() const
    {
        return m_opposingVerticesPerEdge;
    }

protected:
    std::vector<Vector3> m_vertices;
    std::vector<Vector3> m_vertices0;
    std::vector<Int3> m_faces;
    std::vector<Int2> m_edges;
    std::vector<std::vector<OpposingVertex>> m_opposingVerticesPerEdge;
};

}  // namespace spg