#using JuMP
using Ipopt
include("mds_metric.jl")

contains_inf(array) = mapreduce(x->(isinf(x)), &, array)

# A data structure representing the problem data and constraints
type MDSProblem
    num_points::Int
    dist::Array{Float64, 2}
    sphere_radius::Float64
    lads::Array{Int, 1}
    nads::Array{Int, 1}
    nucleolus_radius::Float64

    MDSProblem(n::Int, dist::Array{Float64, 2}) = 
        new(n, dist, 1.0, Array(Int,0), Array(Int,0), 0)
    MDSProblem(n::Int, dist::Array{Float64, 2}, r::Float64) = 
        new(n, dist, r, Array(Int,0), Array(Int,0), 0)
    # New problem with just a radius constraint
    MDSProblem(dist::Array{Float64, 2}, r::Float64) = 
        new(size(dist)[1], dist, r, Array(Int,0), Array(Int,0), 0)
    MDSProblem(dist::Array{Float64, 2}, r::Float64, lads_file::String) = 
        new(size(dist)[1], dist, r, load_lads(lads_file), Array(Int,0), 0)
    MDSProblem(dist::Array{Float64, 2}, r::Float64, lads_file::String, nads_file::String, nucleolus_r::Float64) = 
        new(size(dist)[1], dist, r, load_lads(lads_file), load_lads(nads_file), nucleolus_r)
end


function load_data(filename::String; exponent=1/3, scale=1)
    # loads data, convert to NxN array of "wish distances"
    contact_data = readcsv(filename, Float64)
    #dist = contact_freq_to_distance_matrix(filename)
    dist = scale./(contact_data).^exponent
    return dist
end

function load_lads(filename::String)
    # should return a Vector{Float64}
    # a "lads file" is just the length followed by a list of indices.
    f = open(filename, "r")
    length = int(strip(readline(f)))
    lads = Array(Int, length)
    for (i, line) in enumerate(readlines(f))
        lads[i] = int(strip(line))
    end
    return lads
end

function remove_zeros(dist::Array{Float64, 2})
    # Removes all-zero and all-inf rows/columns from a contact map
    dist_index = 1
    while dist_index <= size(dist)[1]
        if mapreduce(x->(isinf(x) || x==0), &, dist[dist_index, 1:end-1])
            dist = [dist[1:dist_index-1, :], dist[dist_index+1:end, :]]
            dist = [dist[:, 1:dist_index-1] dist[:,dist_index+1:end]]
        else
            dist_index += 1
        end
    end
    return dist
end

function remove_infs(prob::MDSProblem)
    # Removes all rows from the distance matrix that have all infinities.
    dist_index = 1
    dist = prob.dist
    while dist_index <= size(dist)[1]
        if mapreduce(x->(isinf(x) || x==0), &, dist[dist_index, 1:end-1])
            dist = [dist[1:dist_index-1, :], dist[dist_index+1:end, :]]
            dist = [dist[:, 1:dist_index-1] dist[:,dist_index+1:end]]
            new_lads = filter(x->x!=dist_index, prob.lads)
            for (i, lad) in enumerate(new_lads)
                if lad > dist_index
                    new_lads[i] -=1
                end
            end
            prob.lads = new_lads
            new_nads = filter(x->x!=dist_index, prob.nads)
            for (i, nad) in enumerate(new_nads)
                if nad > dist_index
                    new_nads[i] -=1
                end
            end
            prob.nads = new_nads
        else
            dist_index += 1
        end
    end
    prob.dist = dist
    prob.num_points = size(dist)[1]
    return prob
end


function make_ipopt_problem(mds::MDSProblem; radius_constraint=false, lad_constraint=false, nad_constraint=false)
    # Creates an Ipopt problem for the given problem instance

    # 3 vars for eaxh point- x, y, z
    num_vars = 3*mds.num_points
    num_constraints = 0
    n, n2 = size(mds.dist)
    dist = mds.dist
    rad = mds.sphere_radius
    nuc_rad = mds.nucleolus_radius
    initial_x = [rand()*2*rad - rad for i in 1:num_vars]
    # m is the number of constraints
    m = 0
    if radius_constraint
        m = mds.num_points
    end
    #println(m)
    # constraint lower bound
    g_L = zeros(m)
    # constraint upper bound
    g_U = fill(rad^2, m)
    if lad_constraint
        for l in mds.lads
            g_L[l] = 0.9*rad^2
        end
    end
    # TODO: add constraint for nucleolus center
    # nuc_c is the index fo the x-element of the nucleolus coordinate.
    nuc_c = 3*mds.num_points + 1
    num_jac_elements=m*3
    if nad_constraint
        num_vars += 3
        m += 3 + length(mds.nads)
        append!(initial_x, [rand()*2*rad - rad for i in 1:3])
        append!(g_L, zeros(3 + length(mds.nads)))
        append!(g_U, fill((rad - nuc_rad)^2, 3))
        append!(g_U, fill(nuc_rad^2, length(mds.nads)))
        num_jac_elements += 3 + 6*length(mds.nads)
    end

    interactions = Array((Int, Int), 0)
    for i in 1:n-1
        for j in i+1:n
            if dist[i,j] != Inf
                push!(interactions, (i, j))
            end
        end
    end
    println("Number of interactions: ", length(interactions))
    println("Number of variables: ", num_vars)
    println("Number of constraints: ", m)
    println("Number of nonzero Jacobian elements: ", num_jac_elements)

    ###### Callback Functions

    function eval_f(x::Vector{Float64})
        val = 0
        for (i, j) in interactions
            i_c = 3*(i-1) + 1
            j_c = 3*(j-1) + 1
            s = (x[i_c] - x[j_c])^2 + (x[i_c+1] - x[j_c+1])^2 + (x[i_c+2] - x[j_c+2])^2
            w = 1/dist[i,j]^2
            val += w*(sqrt(s) - dist[i,j])^2
        end
        #println("x max: ", maximum(x))
        return val
    end

    function eval_g(x::Vector{Float64}, g::Vector{Float64})
        g[:] = zeros(size(g))
        if !nad_constraint
            for i in 1:m
                i_c = 3*(i-1) + 1
                g[i] = x[i_c]^2 + x[i_c+1]^2 + x[i_c+2]^2
            end
        else
            for i in 1:mds.num_points
                i_c = 3*(i-1) + 1
                g[i] = x[i_c]^2 + x[i_c+1]^2 + x[i_c+2]^2
            end
            g_index = mds.num_points+1
            g[g_index] = x[nuc_c]^2 + x[nuc_c+1]^2+x[nuc_c+2]^2
            g_index += 1
            for i in mds.nads
                i_c = 3*(i-1) + 1
                g[g_index] = (x[i_c] - x[nuc_c])^2 + (x[i_c+1] - x[nuc_c+1])^2
                    + (x[i_c+2] - x[nuc_c+2])^2
                g_index += 1
            end
        end
    end

    function eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64})
        grad_f[:] = zeros(size(grad_f))
        for (i, j) in interactions
            # PROBLEm: we get huge gradients
            i_c = 3*(i-1) + 1
            j_c = 3*(j-1) + 1
            #s = sum((x[i_c:i_c+2] .- x[j_c:j_c+2]).^2)
            s = (x[i_c] - x[j_c])^2 + (x[i_c+1] - x[j_c+1])^2 + (x[i_c+2] - x[j_c+2])^2
            s2 = sqrt(s)
            f = s2-dist[i,j]
            if s2 > 0.0000001
                w = 1/dist[i,j]^2
                common_term = w*(1-dist[i,j]/s2)
                #if abs(common_term) > 50
                #    krintln("index: ", i, j)
                #    println("x: ", x[i_c:i_c+2], " ", x[j_c:j_c+2])
                #end
                grad_f[i_c] += common_term*2*(x[i_c] - x[j_c])
                grad_f[i_c+1] += common_term*2*(x[i_c+1] - x[j_c+1])
                grad_f[i_c+2] += common_term*2*(x[i_c+2] - x[j_c+2])
                grad_f[j_c] -= common_term*2*(x[i_c] - x[j_c])
                grad_f[j_c+1] -= common_term*2*(x[i_c+1] - x[j_c+1])
                grad_f[j_c+2] -= common_term*2*(x[i_c+2] - x[j_c+2])
            end
        end
        #println(grad_f)
        #println("grad max: ", maximum(grad_f))
        #println("grad min: ", minimum(grad_f))
        #mi = indmax(grad_f)
        #println("index of grad max: ", mi)
        #println("x max: ", maximum(x))
        #println("x min: ", minimum(x))
    end

    function eval_jac_g(
        x::Vector{Float64},         # Current solution
        mode,                       # Either :Structure or :Values
        rows::Vector{Int32},        # Sparsity structure - row indices
        cols::Vector{Int32},        # Sparsity structure - column indices
        values::Vector{Float64})
        #println("eval jac g")
        if mode == :Structure
            #println("eval jac g structure")
            if !nad_constraint
                for i in 1:m
                    i_c = 3*(i-1) + 1
                    rows[i_c] = i
                    rows[i_c+1] = i
                    rows[i_c+2] = i
                    cols[i_c] = i_c
                    cols[i_c+1] = i_c+1
                    cols[i_c+2] = i_c+2
                end
            else
                for i in 1:mds.num_points
                    i_c = 3*(i-1) + 1
                    #println(i_c)
                    rows[i_c] = i
                    rows[i_c+1] = i
                    rows[i_c+2] = i
                    cols[i_c] = i_c
                    cols[i_c+1] = i_c+1
                    cols[i_c+2] = i_c+2
                end
                constraint_index = mds.num_points + 1
                i_c = 3*(mds.num_points)+1
                rows[i_c] = constraint_index
                rows[i_c+1] = constraint_index
                rows[i_c+2] = constraint_index
                cols[i_c] = nuc_c
                cols[i_c+1] = nuc_c+1
                cols[i_c+2] = nuc_c+2
                constraint_index += 1
                i_c += 3
                for i in mds.nads
                    #println(i_c)
                    x_index = 3*(i-1)+1
                    rows[i_c] = constraint_index
                    rows[i_c+1] = constraint_index
                    rows[i_c+2] = constraint_index
                    rows[i_c+3] = constraint_index
                    rows[i_c+4] = constraint_index
                    rows[i_c+5] = constraint_index
                    cols[i_c] = x_index
                    cols[i_c+1] = x_index+1
                    cols[i_c+2] = x_index+2
                    cols[i_c+3] = nuc_c
                    cols[i_c+4] = nuc_c+1
                    cols[i_c+5] = nuc_c+2
                    i_c += 6
                    constraint_index += 1
                end
            end
        else
            if !nad_constraint
                for i in 1:m
                    i_c = 3*(i-1) + 1
                    values[i_c] = 2*x[i_c]
                    values[i_c+1] = 2*x[i_c+1]
                    values[i_c+2] = 2*x[i_c+2]
                end
            else
                for i in 1:mds.num_points
                    i_c = 3*(i-1) + 1
                    values[i_c] = 2*x[i_c]
                    values[i_c+1] = 2*x[i_c+1]
                    values[i_c+2] = 2*x[i_c+2]
                end
                #constraint_index = mds.num_points+1
                i_c = 3*(mds.num_points)+1
                values[i_c] = 2*x[nuc_c]
                values[i_c+1] = 2*x[nuc_c+1]
                values[i_c+2] = 2*x[nuc_c+1]
                i_c += 3
                for i in mds.nads
                    x_index = 3*(i-1)+1
                    values[i_c] = 2*(x[x_index]-x[nuc_c])
                    values[i_c+1] = 2*(x[x_index+1]-x[nuc_c+1])
                    values[i_c+2] = 2*(x[x_index+2]-x[nuc_c+2])
                    values[i_c+3] = -2*(x[x_index]-x[nuc_c])
                    values[i_c+4] = -2*(x[x_index+1]-x[nuc_c+1])
                    values[i_c+5] = -2*(x[x_index+2]-x[nuc_c+2])
                    i_c += 6
                end
            end
        end
        #println("finished eval jac g")
    end

    function eval_h(
        x::Vector{Float64},         # Current solution
        mode,                       # Either :Structure or :Values
        rows::Vector{Int32},        # Sparsity structure - row indices
        cols::Vector{Int32},        # Sparsity structure - column indices
        obj_factor::Float64,        # Lagrangian multiplier for objective
        lambda::Vector{Float64},    # Multipliers for each constraint
        values::Vector{Float64})    # The values of the Hessian

        if mode == :Structure
            # rows[...] = ...
            # ...
            # cols[...] = ...
        else
            # values[...] = ...
        end
    end

    function intermediate(
        alg_mod::Int,
        iter_count::Int,
        obj_value::Float64,
        inf_pr::Float64, inf_du::Float64,
        mu::Float64, d_norm::Float64,
        regularization_size::Float64,
        alpha_du::Float64, alpha_pr::Float64,
        ls_trials::Int)
        # ...
        return true  # Keep going
    end

    problem = createProblem(
        num_vars::Int,                     # Number of variables
        fill(-rad, num_vars)::Vector{Float64},       # Variable lower bounds
        fill(rad, num_vars)::Vector{Float64},       # Variable upper bounds
        m::Int,                     # Number of constraints
        g_L,       # Constraint lower bounds
        g_U,       # Constraint upper bounds
        num_jac_elements,              # Number of non-zeros in Jacobian
        0,             # Number of non-zeros in Hessian
        eval_f,                     # Callback: objective function
        eval_g,                     # Callback: constraint evaluation
        eval_grad_f,                # Callback: objective function gradient
        eval_jac_g,                 # Callback: Jacobian evaluation
        nothing)           # Callback: Hessian evaluation
    addOption(problem, "hessian_approximation", "limited-memory")
    addOption(problem, "tol", 1e-1)
    addOption(problem, "acceptable_tol", 1e-1)
    addOption(problem, "max_iter", 500)
    problem.x = initial_x
    return problem


end
