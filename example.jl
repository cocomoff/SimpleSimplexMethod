include("./code.jl")
using Main.SimpleSimplex

# problem example
c = [-3; -2; -1; -5; 0; 0; 0]
A = [7 3 4 1 1 0 0;
     2 1 1 5 0 1 0;
     1 4 5 2 0 0 1]
b = [7; 3; 8]

# BFS
println("BFS")
opt_x, obj = searchBFS(c, A, b)
println(obj)
println(opt_x)
println()

# Simplex method
println("Simplex method")
opt_x, obj = simplex_method(c, A, b)
println(obj)
println(opt_x)