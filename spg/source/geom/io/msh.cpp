#include <spg/geom/io/msh.h>

// Since there is no real support for external obj creations for now, hide this code behind a compilation flag
#ifdef MSHIO_ACTIVE
#include <mshio/mshio.h>
#endif

namespace spg::io
{
std::tuple<std::vector<Vector3>, std::vector<Int4>> loadMsh(const std::string &filePath)
{
    std::vector<Vector3> vertices;
    std::vector<Int4> tets;
#ifdef MSHIO_ACTIVE
    // TODO: Make a more secure load, checking that sizes are the expected ones
    mshio::MshSpec spec = mshio::load_msh(filePath);

    vertices.reserve(spec.nodes.num_nodes);
    for (size_t i = 0; i < spec.nodes.num_entity_blocks; i++) {
        mshio::NodeBlock &block = spec.nodes.entity_blocks[i];
        for (size_t j = 0; j < block.num_nodes_in_block; j++) {
            vertices.emplace_back(block.data[j * 3 + 0], block.data[j * 3 + 1], block.data[j * 3 + 2]);
        }
    }
    std::cout << "nverts " << vertices.size();
    tets.reserve(spec.elements.num_elements);
    for (size_t i = 0; i < spec.elements.num_entity_blocks; i++) {
        mshio::ElementBlock &block = spec.elements.entity_blocks[i];
        const size_t n = mshio::nodes_per_element(block.element_type);
        for (size_t j = 0; j < block.num_elements_in_block; j++) {
            tets.push_back({static_cast<int>(block.data[j * 5 + 1]) - 1,
                            static_cast<int>(block.data[j * 5 + 2]) - 1,
                            static_cast<int>(block.data[j * 5 + 3]) - 1,
                            static_cast<int>(block.data[j * 5 + 4]) - 1});
        }
    }
    std::cout << "nverts " << vertices.size();
#else
    throw std::runtime_error("loadMsh requires compiling with MshIO activated");
#endif
    return {std::move(vertices), std::move(tets)};
}
}  // namespace spg::io