#pragma once

#include <spg/types.h>
#include <tuple>
namespace spg::io
{
std::tuple<std::vector<Vector3>, std::vector<Int3>> loadObj(const std::string &filePath);
}  // namespace spg::io