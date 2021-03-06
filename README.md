# Simplex Julia

- Julia implementation of [Simplex method](https://en.wikipedia.org/wiki/Simplex_algorithm).
- The code is based on the book: [Julia Programming for Operations Research, A Primer on Computing](https://www.softcover.io/read/7b8eb7d0/juliabook), but the julia version is assumed to be 1.5.


# Required packages

- Combinatorics v1.0.2


# How to use

- Run *example.jl* or use Julia terminal as follows.

```julia
using Revise
includet("code.jl")
using Main.SimplexMethod
c = [-3; -2; -1; -5; 0; 0; 0]
A = [7 3 4 1 1 0 0; 2 1 1 5 0 1 0; 1 4 5 2 0 0 1]
b = [7; 3; 8]
simplex_method(c, A, b)
```

# Output of example.jl

```txt
BFS
-5.051724137931034
[0.1724137931034483, 1.8793103448275859, 0.0, 0.15517241379310343, 0.0, 0.0, 0.0]

Simplex method
------+-------------------------------------------------+-------
      |  3.00   2.00   1.00   5.00   0.00   0.00   0.00 |   0.00
------+-------------------------------------------------+-------
x[ 5] |  7.00   3.00   4.00   1.00   1.00   0.00   0.00 |   7.00
x[ 6] |  2.00   1.00   1.00   5.00   0.00   1.00   0.00 |   3.00
x[ 7] |  1.00   4.00   5.00   2.00   0.00   0.00   1.00 |   8.00
------+-------------------------------------------------+-------
Pivoting: entering = x_1, exiting = x_5
------+-------------------------------------------------+-------
      |  0.00   0.71  -0.71   4.57  -0.43   0.00   0.00 |  -3.00
------+-------------------------------------------------+-------
x[ 1] |  1.00   0.43   0.57   0.14   0.14   0.00   0.00 |   1.00
x[ 6] |  0.00   0.14  -0.14   4.71  -0.29   1.00   0.00 |   1.00
x[ 7] |  0.00   3.57   4.43   1.86  -0.14   0.00   1.00 |   7.00
------+-------------------------------------------------+-------
Pivoting: entering = x_2, exiting = x_7
------+-------------------------------------------------+-------
      |  0.00   0.00  -1.60   4.20  -0.40   0.00  -0.20 |  -4.40
------+-------------------------------------------------+-------
x[ 1] |  1.00   0.00   0.04  -0.08   0.16   0.00  -0.12 |   0.16
x[ 6] |  0.00   0.00  -0.32   4.64  -0.28   1.00  -0.04 |   0.72
x[ 2] |  0.00   1.00   1.24   0.52  -0.04   0.00   0.28 |   1.96
------+-------------------------------------------------+-------
Pivoting: entering = x_4, exiting = x_6
------+-------------------------------------------------+-------
      |  0.00   0.00  -1.31   0.00  -0.15  -0.91  -0.16 |  -5.05
------+-------------------------------------------------+-------
x[ 1] |  1.00   0.00   0.03   0.00   0.16   0.02  -0.12 |   0.17
x[ 4] |  0.00   0.00  -0.07   1.00  -0.06   0.22  -0.01 |   0.16
x[ 2] |  0.00   1.00   1.28   0.00  -0.01  -0.11   0.28 |   1.88
------+-------------------------------------------------+-------
-5.051724137931035
[0.17241379310344832, 1.879310344827586, 0.0, 0.15517241379310345, 0.0, 0.0, 0.0]
```