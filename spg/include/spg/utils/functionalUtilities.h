#pragma once

namespace spg
{
// TODO: Review if more std::forward's are needed for better performance
template <typename Func, typename Tuple>
void apply_each(Func &&function, Tuple &&tuple)
{
    std::apply(
        [&func = std::forward<Func>(function)](auto &&...args) { (func(std::forward<decltype(args)>(args)), ...); },
        std::forward<Tuple>(tuple));
}
}  // namespace spg