# WaveFront.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wtegtow.github.io/WaveFront.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wtegtow.github.io/WaveFront.jl/dev/)
[![Build Status](https://github.com/wtegtow/WaveFront.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wtegtow/WaveFront.jl/actions/workflows/CI.yml?query=branch%3Amain)

**WaveFront.jl** solves the **Eikonal equation** of the form:

$$
|\nabla T(x)| = \frac{1}{v(x)}
$$

using the **Fast Sweeping Method**.


### Installation and Quick Start

```julia
using Pkg 
Pkg.add(url="https://github.com/wtegtow/WaveFront.jl")

using WaveFront

# 2D example 
x_coords = 0:20:1000 
y_coords = 0:20:1000
source_coords = [(200, 500), (500, 200)]
velocity = ones(length(x_coords), 
                length(y_coords)) .* 500;

grid = Grid2D(x_coords, y_coords, velocity)
convergence_criterium = 1e-6 # L1 norm for 2D
tt = fast_sweeping(grid, [source_coords], ϵ = 1e-6, max_iter=50, verbose=false);

# 3D example 
z_coords = 0:20:1000 
source_coords = [(500, 500, 500)]
velocity = ones(length(x_coords), 
                length(y_coords),
                length(z_coords)) .* 500;

grid = Grid3D(x_coords, y_coords, z_coords, velocity)
convergence_criterium = 0.1 # Supremum norm for 3D
tt = fast_sweeping(grid, [source_coords], ϵ = convergence_criterium, max_iter=50, verbose=false);

```

Additional examples are provided in the `examples/` folder.