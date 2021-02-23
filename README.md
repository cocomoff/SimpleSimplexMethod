# Simplex Julia

- Julia implementation of [Simplex method](https://en.wikipedia.org/wiki/Simplex_algorithm).
- The code is based on the book: [Julia Programming for Operations Research, A Primer on Computing](https://www.softcover.io/read/7b8eb7d0/juliabook), but the julia version is assumed to be 1.5.


# Required packages

- Combinatorics v1.0.2


# How to use

- Run *example.jl* or use Julia terminal as follows.

```julia
julia> using Revise
julia> includet("code.jl")
julia> using Main.SimplexMethod
julia> c = [-3; -2; -1; -5; 0; 0; 0]
julia> A = [7 3 4 1 1 0 0; 2 1 1 5 0 1 0; 1 4 5 2 0 0 1]
julia> b = [7; 3; 8]
julia> simplex_method(c, A, b)
```