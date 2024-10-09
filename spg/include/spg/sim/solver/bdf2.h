#pragma once

#include <spg/sim/solver/implicitEulerBase.h>

#include <vector>
namespace spg
{
namespace solver
{
// refs: "Stable but Responsive Cloth", 2002 and "Dynamic Deformables: Implementation and Production Practicalities",
// 2022
class BDF2 : public ImplicitEulerBase
{
public:
    virtual void step();
    virtual void reset();

protected:
    VectorX m_xMinusOne;
    VectorX m_vMinusOne;
};
}  // namespace solver
}  // namespace spg