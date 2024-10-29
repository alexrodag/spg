#pragma once

#include <list>
#include <vector>

namespace spg
{
namespace coloring
{
// Very simple graph class and greedy graph coloring strategy
// Ref: https://en.wikipedia.org/wiki/Greedy_coloring
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

// Color stencils by adjacency through vertices
std::vector<int> colorStencils(const FlatStencils &flatStencils);

// Color vertices by adjacency through stencils
std::vector<int> colorVertices(int numVertices, const std::vector<FlatStencils> &flatStencilsSet);
}  // namespace coloring
}  // namespace spg