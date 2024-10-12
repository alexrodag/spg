# SPG

SPG (short for Spring) is a research-focused C++ physics simulation library for easy implementation and comparison of solvers and mechanical energy models. It provides a common interface for solvers and energies, making all solver implementations compatible with all energies without requiring additional effort. A set of solvers and energies are already implemented in SPG. I personally use it to implement and compare interesting ideas from state-of-the-art research.

![SPG solver comparison example](media/spg-solver-compare.gif)

General nonlinear energies $E(x) = f(x)$ are supported. Energies of the particular form $E(x) = \frac{1}{2}C(x)^T K C(x)$, with $C(x)$ a vector of nonlinear constraint functions and $K$ a stiffness matrix, can also be defined by providing just the implementation of $C(x)$ and become automatically compatible with constraint-based solvers (e.g. XPBD) as well as energy-based and force-based solvers.

Through the use of autodiff and templates (using [TinyAD](https://github.com/patr-schm/TinyAD) as autodiff backend), automatic first and second order derivatives are provided from a simple energy or constraint definition. This helps to keep the code to a minimum, very close to the actual formulas of the corresponding solvers and models. It is also possible to override these automatic versions with more efficient custom implementations if needed.

As an example, by providing the following constraint code for the [Discrete Bending Energy](https://media.disneyanimation.com/uploads/production/publication_asset/75/asset/WDAS_TR_201307.pdf), the corresponding constraint derivatives, energy, forces and force jacobians are automatically available, making it compatible with all the implemented solvers:

``` C++
    [...]
    const Vector3T<ADouble> e0 = x1 - x0;
    const Vector3T<ADouble> e3 = x2 - x1;
    const Vector3T<ADouble> e4 = x3 - x1;

    const Vector3T<ADouble> n1 = e0.cross(e3).normalized();
    const Vector3T<ADouble> n2 = -e0.cross(e4).normalized();

    const auto theta = atan2((n1.cross(n2)).dot(e0.normalized()), n1.dot(n2));
    const auto constraint = theta - restTheta;
    [...]
```

### Current implemented solvers
- [XPBD](https://matthias-research.github.io/pages/publications/XPBD.pdf)
- [VBD (Vertex Block Descent)](https://arxiv.org/pdf/2403.06321) - [Author's source code](https://github.com/AnkaChan/Gaia)
- [Baraff-Witkin Implicit Euler](https://www.cs.cmu.edu/~baraff/papers/sig98.pdf)
- [Robust Newton Backward Euler](https://drive.google.com/file/d/1KbRVF7fk5AonJelIcivFruS2cG3ke-Oi/view)
- [BDF2](https://www.tkim.graphics/DYNAMIC_DEFORMABLES/DynamicDeformables.pdf)
- [Preconditioned Gradient Descent](https://wanghmin.github.io/publication/wang-2016-dme/Wang-2016-DME.pdf)
- [Quasistatic equilibrium](https://pcs-sim.github.io/)
- Static equilibrium

### Current implemented energies
- Spring
- [Strain-based spring](https://cg.informatik.uni-freiburg.de/publications/2007_SCA_ropes.pdf)
- [Baraff-Witkin bending](https://www.cs.cmu.edu/~baraff/papers/sig98.pdf)
- [Quadratic bending](https://cims.nyu.edu/gcl/papers/bergou2006qbm.pdf)
- [Discrete bending](https://media.disneyanimation.com/uploads/production/publication_asset/75/asset/WDAS_TR_201307.pdf)
- [Baraff-Witkin membrane](https://www.cs.cmu.edu/~baraff/papers/sig98.pdf)
- [Choi membrane](https://diglib.eg.org/server/api/core/bitstreams/84fd4836-05cb-4da6-a230-d965b93335a2/content)
- [StVK membrane](https://inria.hal.science/inria-00394466/PDF/tensile.pdf)
- [Stable Neo-Hookean](https://graphics.pixar.com/library/StableElasticity/paper.pdf)
- [StVK](https://en.wikipedia.org/wiki/Hyperelastic_material#Saint_Venant%E2%80%93Kirchhoff_model)

The current implementations are focused on parallel CPU simulation of deformables, including cloth, spring systems and finite element models. Some of the solvers meant for GPU computation may underperform, but this common framework allows to analyse and compare other important aspects such as convergence and stability properties.

This is still a experimental project, which is the best excuse I can think of to justify the lack of proper documentation and unit testing.

I intend to explore several big features in the future, in no particular order (if you are interested in contributing, let me know):
- Unified energy damping model
- Proper boundary conditions
- Collision detection
- Rigid bodies
- Differentiable simulation
- GPU support

## Compilation and dependencies
CMake is used to generate and compile the project. It has been tested in Windows 11 with VSCode and MSVC 19, and in Ubuntu 24.04. You can configure it in your preferred IDE. The `apps/` folder contains demo examples, including a small elastic pendulum simulation with the basics of setting up and simulating a scene, and a solver comparison demo that allows to set up different scenes and simulate them with multiple configurable solvers. 

SPG uses [Eigen](https://eigen.tuxfamily.org/) and [TinyAD](https://github.com/patr-schm/TinyAD) as dependencies. The solver comparison demo also uses [Polyscope](https://polyscope.run/) for GUI and rendering. All of these come as submodules.

For a quick Terminal compilation on both Windows or Ubuntu, you can run:

```
git clone https://github.com/alexrodag/spg.git
cd spg
git submodule update --init --recursive
mkdir build
cd build
cmake ..
cmake --build . --config Release -j
```

**Note**: In Ubuntu, if there is an error when configuring polyscope, you may need to run the following (check updated instructions in https://polyscope.run/building/)
)
```
sudo apt install xorg-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev
```

---

Author: [Alex Rodriguez](https://sites.google.com/view/alejandrora)