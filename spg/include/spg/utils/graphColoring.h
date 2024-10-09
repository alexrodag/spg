#pragma once

#include <list>
#include <vector>

namespace spg
{
namespace coloring
{
// Very simple graph class and graph coloring strategy
class Graph
{
public:
    void setNumNodes(int numNodes);
    void addEdge(int i, int j);

    // Assigns colors (starting from 0) to all nodes
    std::vector<int> greedyColoring();

protected:
    std::vector<std::list<int>> m_adjacency;
};

struct FlatStencils {
    std::vector<int> entries;
    int stencilSize;
};

// Color stencils by adjacency
std::vector<int> colorStencils(const FlatStencils &flatStencils);
std::vector<int> colorVertices(int numVertices, const std::vector<FlatStencils> &flatStencilsSet);
}  // namespace coloring
}  // namespace spg