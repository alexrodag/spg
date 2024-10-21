#pragma once

#include <spg/types.h>
#include <tuple>
namespace spg::io
{
std::tuple<std::vector<Vector3>, std::vector<Int4>> loadMsh(const std::string &filePath);
}  // namespace spg::io