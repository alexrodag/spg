#include <spg/utils/graphColoring.h>
#include <stdexcept>

#include <algorithm>

namespace spg::coloring
{
void Graph::setNumNodes(const int numNodes)
{
    m_adjacency.resize(numNodes);
}
void Graph::addEdge(const int i, const int j)
{
    if (i >= m_adjacency.size() || j >= m_adjacency.size()) {
        throw std::runtime_error("Trying to add edge to graph with invalid node indices");
    }
    m_adjacency[i].push_back(j);
    m_adjacency[j].push_back(i);
}

std::vector<int> Graph::greedyColoring()
{
    const int nbNodes = static_cast<int>(m_adjacency.size());
    std::vector<int> nodeColors(nbNodes, -1);

    // Assign the first color to first node
    nodeColors[0] = 0;

    std::vector<char> colorIsAvailable(nbNodes);

    // Assign colors to remaining nodes
    for (int n = 1; n < nbNodes; n++) {
        std::fill(colorIsAvailable.begin(), colorIsAvailable.end(), 1);
        // Flag unavailable colors
        std::for_each(m_adjacency[n].begin(), m_adjacency[n].end(), [&colorIsAvailable, &nodeColors](const int &n) {
            if (nodeColors[n] != -1) {
                colorIsAvailable[nodeColors[n]] = 0;
            }
        });
        // Assign the first available color
        for (int color = 0; color < nbNodes; ++color) {
            if (colorIsAvailable[color] == 1) {
                nodeColors[n] = color;
                break;
            }
        }
    }
    return nodeColors;
}

std::vector<int> colorStencils(const FlatStencils &flatStencils)
{
    Graph graph;
    std::vector<std::list<int>> stencilsPerVertex;
    for (int i = 0; i < flatStencils.entries.size(); ++i) {
        const int vertexId = flatStencils.entries[i];
        const int stencilId = i / flatStencils.stencilSize;
        stencilsPerVertex.resize(std::max(stencilsPerVertex.size(), static_cast<size_t>(vertexId + 1)));
        stencilsPerVertex[vertexId].push_back(stencilId);
    }
    graph.setNumNodes(static_cast<int>(flatStencils.entries.size()) / flatStencils.stencilSize);
    for (const auto &vertexStencils : stencilsPerVertex) {
        for (auto it1 = vertexStencils.begin(); it1 != vertexStencils.end(); it1++) {
            for (auto it2 = std::next(it1); it2 != vertexStencils.end(); it2++) {
                graph.addEdge(*it1, *it2);
            }
        }
    }
    return graph.greedyColoring();
}

std::vector<int> colorVertices(int numVertices, const std::vector<FlatStencils> &flatStencilsSet)
{
    Graph graph;
    graph.setNumNodes(numVertices);
    for (auto &flatStencils : flatStencilsSet) {
        for (int i = 0; i < flatStencils.entries.size(); ++i) {
            for (int j = i + 1; j % flatStencils.stencilSize != 0; ++j) {
                graph.addEdge(flatStencils.entries[i], flatStencils.entries[j]);
            }
        }
    }
    return graph.greedyColoring();
}
}  // namespace spg::coloring