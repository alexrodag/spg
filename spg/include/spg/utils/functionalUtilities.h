#pragma once

namespace spg
{
// TODO: Review if more std::forward's are needed for better performance
template <typename Func, typename Tuple>
void apply_each(Func &&function, Tuple &&tuple)
{
    std::apply([&function](auto &&...args) { (function(args), ...); }, std::forward<Tuple>(tuple));
}
}  // namespace spg