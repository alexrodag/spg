#pragma once

#include <spg/types.h>
#include <set>
#include <stdexcept>

namespace spg
{
class TetrahedralMesh
{
public:
    TetrahedralMesh(const std::vector<Vector3> &vertices,
                    const std::vector<Vector3> &vertices0,
                    const std::vector<Int4> &tets)
        : m_vertices(vertices)
        , m_vertices0(vertices0)
        , m_tets(tets)
    {
        const int nVertices{static_cast<int>(m_vertices.size())};
        if (std::find_if(m_tets.begin(), m_tets.end(), [nVertices](const auto &tet) {
                return std::any_of(
                    tet.begin(), tet.end(), [nVertices](const auto &index) { return index < 0 || index >= nVertices; });
            }) != m_tets.end()) {
            throw std::runtime_error("Tetrahedral mesh cell contains index to invalid vertex");
        }
        if (std::find_if(m_tets.begin(), m_tets.end(), [nVertices](const auto &tet) {
                return tet[0] == tet[1] || tet[0] == tet[2] || tet[0] == tet[3] || tet[1] == tet[2] ||
                       tet[1] == tet[3] || tet[2] == tet[3];
            }) != m_tets.end()) {
            throw std::runtime_error("Tetrahedral mesh cell contains repeated vertex index");
        }
        // Compute edges info
        std::set<Int2> edges;
        for (int t = 0; t < m_tets.size(); ++t) {
            const auto &tet{m_tets[t]};
            for (int i = 0; i < 3; ++i) {
                for (int j = i + 1; j < 4; ++j) {
                    const int id0 = tet[i];
                    const int id1 = tet[j];
                    edges.insert({std::min(id0, id1), std::max(id0, id1)});
                }
            }
        }
        for (const auto &edge : edges) {
            m_edges.emplace_back(edge);
        }
    }
    const std::vector<Vector3> &vertices() const { return m_vertices; }
    const std::vector<Vector3> &vertices0() const { return m_vertices0; }
    const std::vector<Int4> &tets() const { return m_tets; }
    const std::vector<Int2> &edges() const { return m_edges; }

protected:
    std::vector<Vector3> m_vertices;
    std::vector<Vector3> m_vertices0;
    std::vector<Int4> m_tets;
    std::vector<Int2> m_edges;
};

}  // namespace spg