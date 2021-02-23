module SimpleSimplex

using LinearAlgebra
using Combinatorics
using Printf

export searchBFS, simplex_method


isnonnegative(x) = length(x[x .< 0]) == 0


function searchBFS(c, A, b; debug_print=false)
    """
    Search all feasible basis using Combinatorics.jl
    """
    obj = Inf
    m, n = size(A)
    @assert rank(A) == m
    opt_x = zeros(n)

    for b_idx in combinations(1:n, m)
        B = A[:, b_idx]
        c_B = c[b_idx]
        x_B = inv(B) * b

        if isnonnegative(x_B)
            z = dot(c_B, x_B)
            if z < obj
                obj = z
                opt_x = zeros(n)
                opt_x[b_idx] = x_B
            end
        end

        if debug_print
            println("Basis:", b_idx)
            println("\t x_B = ", x_B)
            println("\t nonnegative? ", isnonnegative(x_B))
            if isnonnegative(x_B)
                println("\t obj = ", dot(c_B, x_B))
            end
        end
    end
    return opt_x, obj
end

mutable struct SimplexTableau
    z_c::Array{Float64}
    Y::Array{Float64}
    x_B::Array{Float64}
    obj::Float64
    b_idx::Array{Int64}
end

function initial_BFS(A, b)
    """
    Find an initial feasible basis for (A, b)
    """
    m, n = size(A)
    comb = collect(combinations(1:n, m))
    for i in length(comb):-1:1
        b_idx = comb[i]
        B = A[:, b_idx]
        x_B = inv(B) * b
        isnonnegative(x_B) && return b_idx, x_B, B
    end
    error("infeasible")
end


function initialize(c, A, b)
    """
    Build an initial tableau
    """
    c = Array{Float64}(c)
    A = Array{Float64}(A)
    b = Array{Float64}(b) 
    m, n = size(A)
  
    # 初期基底
    b_idx, x_B, B = initial_BFS(A,b)
    Y = inv(B) * A
    c_B = c[b_idx]
    obj = dot(c_B, x_B)
  
    z_c = zeros(1, n)
    n_idx = setdiff(1:n, b_idx)
    z_c[n_idx] = c_B' * inv(B) * A[:,n_idx] - c[n_idx]'
  
    SimplexTableau(z_c, Y, x_B, obj, b_idx)
end


function print_tableau(t::SimplexTableau)
    """
    Print the given tableau t to terminal
    """
    m, n = size(t.Y)
    hline0 = repeat("-", 6)
    hline1 = repeat("-", 7*n)
    hline2 = repeat("-", 7)
    hline = join([hline0, "+", hline1, "+", hline2])

    println(hline)

    Printf.@printf("%6s|", "")
    for j in 1:length(t.z_c)
        Printf.@printf("%6.2f ", t.z_c[j])
    end
    Printf.@printf("| %6.2f\n", t.obj)

    Printf.println(hline)

    for i in 1:m
        Printf.@printf("x[%2d] |", t.b_idx[i])
        for j in 1:n
            Printf.@printf("%6.2f ", t.Y[i,j])
        end
        Printf.@printf("| %6.2f\n", t.x_B[i])
    end
    println(hline)
end


isOptimal(tableau) = findfirst( tableau.z_c .> 0 ) === nothing


function pivot_point(t::SimplexTableau)
    """
    Finding a pivot (entering/exiting indices)
    """

    # entering index
    entering = findfirst(t.z_c .> 0)
    (entering == nothing) && error("Optimal")
    entering = LinearIndices(t.z_c)[entering]

    # exiting index
    pos_idx = findall( t.Y[:, entering] .> 0 )
    (length(pos_idx) == 0) && error("Unbounded")
    exiting = pos_idx[ argmin( t.x_B[pos_idx] ./ t.Y[pos_idx, entering])]
    return entering, exiting
end


function pivoting!(t::SimplexTableau)
    """
    Pivot operation on the given tableau t
    """
    m, n = size(t.Y)

    entering, exiting = pivot_point(t)
    println("Pivoting: entering = x_$entering, exiting = x_$(t.b_idx[exiting])")

    # Pivoting: exiting-row, entering-column
    # updating exiting-row
    coef = t.Y[exiting, entering]
    t.Y[exiting, :] /= coef
    t.x_B[exiting] /= coef

    # updating other rows of Y
    for i in setdiff(1:m, exiting)
        coef = t.Y[i, entering]
        t.Y[i, :] -= coef * t.Y[exiting, :]
        t.x_B[i] -= coef * t.x_B[exiting]
    end

    # updating the row for the reduced costs
    coef = t.z_c[entering]
    t.z_c -= coef * t.Y[exiting, :]'
    t.obj -= coef * t.x_B[exiting]

    # Updating b_idx
    t.b_idx[ findall(t.b_idx.==t.b_idx[exiting]) ] .= entering
end

function simplex_method(c, A, b; debug_print=false)
    """
    Simplex method
    """
    tableau = initialize(c, A, b)
    print_tableau(tableau)
    
    while !isOptimal(tableau)
        pivoting!(tableau)
        print_tableau(tableau)
    end
    
    opt_x = zeros(length(c))
    opt_x[tableau.b_idx] = tableau.x_B
    opt_x, tableau.obj
end


end