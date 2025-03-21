#include <spg/geom/io/obj.h>

// Since there is no real support for external obj creations for now, hide this code behind a compilation flag
#ifdef TINYOBJLOADER_ACTIVE
#include <tiny_obj_loader.h>
#endif

namespace spg::io
{
std::tuple<std::vector<Vector3>, std::vector<Int3>> loadObj(const std::string &filePath)
{
#ifdef TINYOBJLOADER_ACTIVE
    std::vector<Vector3> vertices;
    std::vector<Int3> faces;
    // TODO: Make a more secure load, checking submeshes, per-face sizes and so on
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::string warn;
    std::string err;
    tinyobj::LoadObj(&attrib, &shapes, nullptr, &warn, &err, filePath.c_str());
    vertices.reserve(attrib.vertices.size() / 3);
    for (int i = 0; i < attrib.vertices.size() / 3; ++i) {
        vertices.push_back({attrib.vertices[3 * i + 0], attrib.vertices[3 * i + 1], attrib.vertices[3 * i + 2]});
    }
    faces.reserve(shapes.front().mesh.indices.size() / 3);
    for (int i = 0; i < shapes.front().mesh.indices.size() / 3; ++i) {
        faces.push_back({shapes.front().mesh.indices[3 * i + 0].vertex_index,
                         shapes.front().mesh.indices[3 * i + 1].vertex_index,
                         shapes.front().mesh.indices[3 * i + 2].vertex_index});
    }
    return {std::move(vertices), std::move(faces)};
#else
    throw std::runtime_error("loadObj requires compiling with tinyobjloader activated, filepath " + filePath +
                             " unused");
#endif
}
}  // namespace spg::io