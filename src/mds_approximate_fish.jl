contains_inf(array) = mapreduce(x->(isinf(x)), &, array)

# A data structure representing the problem data and constraints
type MDSProblem
    num_points::Int
    dist::Array{Float64, 2}
    sphere_radius::Float64
    # constraint and control loci are locations of fish probes,
    # constraint_ratios are ratios of locus distances to the "base" distance
    control_locus_1::Int
    control_locus_2::Int
    c_loci_1::Array{Int, 1}
    c_loci_2::Array{Int, 1}
    constraint_ratios::Array{Float64,1}
    # original distance between control_locus_1 and control_locus_2
    c1_c2_distance::Float64

    # Constructors
    # New problem with just a distance matrix
    MDSProblem(dist::Array{Float64, 2}) = 
        new(size(dist)[1], dist, 1.0, 0, 0, 
            Array(Int,0), Array(Int,0), Array(Float64,0), 0.0)
    # New problem with just a radius constraint
    MDSProblem(dist::Array{Float64, 2}, r::Float64) = 
        new(size(dist)[1], dist, r, 0, 0, 
            Array(Int,0), Array(Int,0), Array(Float64,0), 0.0)
    # New problem with all constraints
    MDSProblem(dist::Array{Float64, 2}, r::Float64,
               cl1::Int,
               cl2::Int,
               c_loci_1::Array{Int,1},
               c_loci_2::Array{Int,1},
               c_ratios::Array{Float64,1},
               d1::Float64) = 
        new(size(dist)[1], dist, r, cl1, cl2, c_loci_1, c_loci_2, c_ratios, d1)
end

function new_mds(dist::Array{Float64,2}, r::Float64, constraints_file::AbstractString;
            resolution::Int=200000, infer_scale::Bool=false)
    # returns a MDSProblem with the constraints in constraints_file
    c = load_constraints(constraints_file, resolution=resolution)
    return MDSProblem(dist, r, c[1], c[2], c[3], c[4], c[5], c[6])
end

function infer_scale(mds::MDSProblem)
    # Infers the scaling factor for a constrained problem, 
    # assuming that the original scaling factor was 1 for the distances.
    d1 = mds.c1_c2_distance
    return d1/mds.dist[mds.control_locus_1, mds.control_locus_2]
end


function load_data(filename::AbstractString; exponent=1/3, scale=1)
    # loads data (in the format of a contact map),
    # convert to NxN array of "wish distances"
    contact_data = readcsv(filename, Float64)
    dist = scale./(contact_data).^exponent
    return dist
end

function load_constraints(filename::AbstractString; resolution::Int=200000)
    # returns a tuple containing all the constraint information - 
    # control_locus_1, control_locus_2, c_loci_1, c_loci_2,
    # constraint_ratios
    #
    # file format:
    # num_constraints
    # locus1    locus2  distance
    # locus1    locus3  distance
    # etc.
    # locus1, locus2, etc. are all bin numbers
    f = open(filename, "r")
    lines = readlines(f)
    n = parse(Int, strip(lines[1])) - 1
    l1 = split(strip(lines[2]))
    c1 = parse(Int, l1[1])
    c2 = parse(Int, l1[2])
    d1 = float(l1[3])
    cl1 = Array(Int, n)
    cl2 = Array(Int, n)
    cr = Array(Float64, n)
    for i in 1:n
        line = lines[2+i]
        l = split(strip(line))
        cl1[i] = int(l[1])
        cl2[i] = int(l[2])
        cr[i] = float(l[3])/d1
    end
    return (c1, c2, cl1, cl2, cr, d1)
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

function get_removed_indices(prob::MDSProblem)
    # Returns a list containing all the indices that are all infs
    dist_index = 1
    dist = prob.dist
    all_removed_indices = Array(Int, 0)
    while dist_index <= size(dist, 1)
        if mapreduce(x->(isinf(x) || x==0), &, dist[dist_index, 1:end-1])
            append!(all_removed_indices, [dist_index])
        end
        dist_index += 1
    end
    return all_removed_indices
end

function print_constraints(prob::MDSProblem)
    # Prints the problem constraints
    println("Constraints")
    println("Locus 1: ", prob.control_locus_1)
    println("Locus 2: ", prob.control_locus_2)
    println("Other loci: ", [(x1, x2) for (x1, x2) in 
            zip(prob.c_loci_1, prob.c_loci_2)])
    println("Length ratios: ", prob.constraint_ratios)
end

function remove_infs(prob::MDSProblem)
    # Removes all rows from the distance matrix that have all infinities.
    # Also shifts the constraints down.
    dist_index = 1
    dist = prob.dist
    while dist_index <= size(dist)[1]
        if mapreduce(x->(isinf(x) || x==0), &, dist[dist_index, 1:end-1])
            dist = [dist[1:dist_index-1, :]; dist[dist_index+1:end, :]]
            dist = [dist[:, 1:dist_index-1] dist[:,dist_index+1:end]]
            # update fish constraints
            if prob.control_locus_1 != 0
                if prob.control_locus_1 >= dist_index
                    prob.control_locus_1 -=1
                end
                if prob.control_locus_2 >= dist_index
                    prob.control_locus_2 -=1
                end
                for i in range(1, length(prob.c_loci_1))
                    if prob.c_loci_1[i] >= dist_index
                        prob.c_loci_1[i] -= 1
                    end
                    if prob.c_loci_2[i] >= dist_index
                        prob.c_loci_2[i] -= 1
                    end
                end
            end
        else
            dist_index += 1
        end
    end
    prob.dist = dist
    prob.num_points = size(dist)[1]
    return prob
end


function make_ipopt_problem(mds::MDSProblem; radius_constraint=false, 
        fish_constraint=false,
        initial_coords=false)
    # Creates an Ipopt problem for the given problem instance

    # 3 vars for each point- x, y, z
    num_vars = 3*mds.num_points
    num_constraints = 0
    n, n2 = size(mds.dist)
    dist = mds.dist
    rad = mds.sphere_radius
    #nuc_rad = mds.nucleolus_radius
    # TODO: add options for initial coordinates
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
    num_jac_elements=m*3
    # the index of g at which fish constraints start
    fish_constraint_index = 0
    # dict of (locus_1, locus_2) => num_constraint
    fish_constraint_dict = Dict{Tuple{Int32, Int32}, Int32}()

    # for FISH constraints, do something
    if fish_constraint
        fish_constraint_index = m+1
        m = m + length(mds.constraint_ratios)
        new_g_L = zeros(length(mds.constraint_ratios))
        new_g_U = zeros(length(mds.constraint_ratios))
        for (i, (i1, i2, r)) in enumerate(zip(mds.c_loci_1, mds.c_loci_2, mds.constraint_ratios))
            # reasoning: if i1, i2 differ from c1, c2, then there are 4 points
            # or 12 variables.
            if haskey(fish_constraint_dict, (i1, i2))
                continue
            end
            n_jac = 12
            if i1 == mds.control_locus_1
                n_jac -= 3
            end
            if i2 == mds.control_locus_1
                n_jac -= 3
            end
            if i1 == mds.control_locus_2
                n_jac -= 3
            end
            if i2 == mds.control_locus_2
                n_jac -= 3
            end
            fish_constraint_dict[i1,i2] = n_jac
            num_jac_elements += n_jac
            # WARNING: arbitrary constant alert 
            new_g_L[i] = -0.00001
            new_g_U[i] = 0.00001
        end
        append!(g_L, new_g_L)
        append!(g_U, new_g_U)
    end


    interactions = Array(Tuple{Int, Int}, 0)
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

    ###### Callback Functions ###########################################

    function eval_f(x::Vector{Float64})
        val = 0
        for (i, j) in interactions
            i_c = 3*(i-1) + 1
            j_c = 3*(j-1) + 1
            s = (x[i_c] - x[j_c])^2 + 
                (x[i_c+1] - x[j_c+1])^2 + 
                (x[i_c+2] - x[j_c+2])^2
            # w is some scaling factor
            w = 1/(dist[i,j])^2
            val += w*(sqrt(s) - dist[i,j])^2
        end
        #println("x max: ", maximum(x))
        return val
    end

    function eval_g(x::Vector{Float64}, g::Vector{Float64})
        # Evaluating constraint functions
        g[:] = zeros(size(g))
        if radius_constraint
            for i in 1:mds.num_points
                i_c = 3*(i-1) + 1
                g[i] = x[i_c]^2 + x[i_c+1]^2 + x[i_c+2]^2
            end
        end
        if fish_constraint
            c1 = 3*(mds.control_locus_1-1)+1
            c2 = 3*(mds.control_locus_2-1)+1
            control_distance = (x[c1]-x[c2])^2 +
                       (x[c1+1]-x[c2+1])^2 + 
                       (x[c1+2]-x[c2+2])^2
            #println("control_distance: ", control_distance)
            for i in fish_constraint_index:m
                c_i = i-fish_constraint_index + 1
                i1 = mds.c_loci_1[c_i] 
                i2 = mds.c_loci_2[c_i] 
                i1_c = 3*(i1-1)+1
                i2_c = 3*(i2-1)+1
                g[i] = (x[i1_c]-x[i2_c])^2 +
                       (x[i1_c+1]-x[i2_c+1])^2 + 
                       (x[i1_c+2]-x[i2_c+2])^2
                g[i] = g[i] - mds.constraint_ratios[c_i]^2*control_distance
                #println(g[i])
            end
            #println(g[fish_constraint_index:m])
        end
    end

    function eval_grad_f(x::Vector{Float64}, grad_f::Vector{Float64})
        # evaluates gradient of the objective function
        grad_f[:] = zeros(size(grad_f))
        for (i, j) in interactions
            i_c = 3*(i-1) + 1
            j_c = 3*(j-1) + 1
            s = (x[i_c] - x[j_c])^2 + (x[i_c+1] - x[j_c+1])^2 + (x[i_c+2] - x[j_c+2])^2
            s2 = sqrt(s)
            f = s2-dist[i,j]
            common_term = 1-dist[i,j]/s2
            w = 1/dist[i,j]^2
            common_term *= w
            grad_f[i_c] += common_term*2*(x[i_c] - x[j_c])
            grad_f[i_c+1] += common_term*2*(x[i_c+1] - x[j_c+1])
            grad_f[i_c+2] += common_term*2*(x[i_c+2] - x[j_c+2])
            grad_f[j_c] -= common_term*2*(x[i_c] - x[j_c])
            grad_f[j_c+1] -= common_term*2*(x[i_c+1] - x[j_c+1])
            grad_f[j_c+2] -= common_term*2*(x[i_c+2] - x[j_c+2])
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
            if fish_constraint
                i_new = fish_constraint_index
                i_c = 3*(i_new-1)+1
                for (i, (i1, i2, r)) in enumerate(zip(mds.c_loci_1, mds.c_loci_2, mds.constraint_ratios))
                    i_c_start = i_c
                    #println("i_c start: ", i_c)
                    n_jac = fish_constraint_dict[i1, i2]
                    cols[i_c] = 3*(i1-1) + 1
                    cols[i_c+1] = 3*(i1-1) + 2
                    cols[i_c+2] = 3*(i1-1) + 3
                    i_c += 3
                    cols[i_c] = 3*(i2-1) + 1
                    cols[i_c+1] = 3*(i2-1) + 2
                    cols[i_c+2] = 3*(i2-1) + 3
                    i_c += 3
                    if i1 != mds.control_locus_1 && i2 != mds.control_locus_1
                        m1 = mds.control_locus_1
                        cols[i_c] = 3*(m1-1) + 1
                        cols[i_c+1] = 3*(m1-1) + 2
                        cols[i_c+2] = 3*(m1-1) + 3
                        i_c += 3
                    end
                    if i1 != mds.control_locus_2 && i2 != mds.control_locus_2
                        m2 = mds.control_locus_2
                        cols[i_c] = 3*(m2-1) + 1
                        cols[i_c+1] = 3*(m2-1) + 2
                        cols[i_c+2] = 3*(m2-1) + 3
                        i_c += 3
                    end
                    #println("i_c: ", i_c)
                    for j in 0:n_jac-1
                        rows[i_c_start+j] = fish_constraint_index+i-1
                    end
                end
            end
            if radius_constraint
                for i in 1:mds.num_points
                    i_c = 3*(i-1) + 1
                    rows[i_c] = i
                    rows[i_c+1] = i
                    rows[i_c+2] = i
                    cols[i_c] = i_c
                    cols[i_c+1] = i_c+1
                    cols[i_c+2] = i_c+2
                end
            end
        else
            if fish_constraint
                i_new = fish_constraint_index
                i_c = 3*(i_new-1)+1
                m1 = mds.control_locus_1
                m1_c = 3*(m1-1)+1
                m2 = mds.control_locus_2
                m2_c = 3*(m2-1)+1
                for (i, (i1, i2, r)) in enumerate(zip(mds.c_loci_1, mds.c_loci_2, mds.constraint_ratios))
                    # 1. gradients for i1
                    i1_c = 3*(i1-1) + 1
                    i2_c = 3*(i2-1) + 1
                    values[i_c] = 2*(x[i1_c]-x[i2_c])
                    values[i_c+1] = 2*(x[i1_c+1]-x[i2_c+1])
                    values[i_c+2] = 2*(x[i1_c+2]-x[i2_c+2])
                    if i1 == mds.control_locus_1
                        values[i_c] -= r^2*2*(x[i1_c]-x[m2_c])
                        values[i_c+1] -= r^2*2*(x[i1_c+1]-x[m2_c+1])
                        values[i_c+2] -= r^2*2*(x[i1_c+2]-x[m2_c+2])
                    elseif i1 == mds.control_locus_2
                        values[i_c] += r^2*2*(x[m1_c] - x[i1_c])
                        values[i_c+1] += r^2*2*(x[m1_c+1] - x[i1_c+1])
                        values[i_c+2] += r^2*2*(x[m1_c+2] - x[i1_c+2])
                    end
                    # 2. gradients for i2
                    i_c += 3
                    values[i_c] = -2*(x[i1_c]-x[i2_c])
                    values[i_c+1] = -2*(x[i1_c+1]-x[i2_c+1])
                    values[i_c+2] = -2*(x[i1_c+2]-x[i2_c+2])
                    if i2 == mds.control_locus_1
                        values[i_c] -= r^2*2*(x[i2_c]-x[m2_c])
                        values[i_c+1] -= r^2*2*(x[i2_c+1]-x[m2_c+1])
                        values[i_c+2] -= r^2*2*(x[i2_c+2]-x[m2_c+2])
                    elseif i2 == mds.control_locus_2
                        values[i_c] += r^2*2*(x[m1_c] - x[i2_c])
                        values[i_c+1] += r^2*2*(x[m1_c+1] - x[i2_c+1])
                        values[i_c+2] += r^2*2*(x[m1_c+2] - x[i2_c+2])
                    end
                    i_c += 3
                    if i1 != mds.control_locus_1 && i2 != mds.control_locus_1
                        # do mds.control_locus_1
                        values[i_c] = -1*r^2*2*(x[m1_c]-x[m2_c])
                        values[i_c+1] = -1*r^2*2*(x[m1_c+1]-x[m2_c+1])
                        values[i_c+2] = -1*r^2*2*(x[m1_c+2]-x[m2_c+2])
                        i_c += 3
                    end
                    if i1 != mds.control_locus_2 && i2 != mds.control_locus_2
                        # do mds.control_locus_2
                        values[i_c] = r^2*2*(x[m1_c] - x[m2_c])
                        values[i_c+1] = r^2*2*(x[m1_c+1] - x[m2_c+1])
                        values[i_c+2] = r^2*2*(x[m1_c+2] - x[m2_c+2])
                        i_c += 3
                    end          
                end
            end
            if radius_constraint
                for i in 1:mds.num_points
                    i_c = 3*(i-1) + 1
                    values[i_c] = 2*x[i_c]
                    values[i_c+1] = 2*x[i_c+1]
                    values[i_c+2] = 2*x[i_c+2]
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
        fill(-Inf, num_vars)::Vector{Float64},       # Variable lower bounds
        fill(Inf, num_vars)::Vector{Float64},       # Variable upper bounds
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
    #addOption(problem, "acceptable_tol", 1e-1)
    addOption(problem, "max_iter", 500)
    addOption(problem, "constr_viol_tol", 1e10)
    addOption(problem, "mu_init", 0.0001)
    addOption(problem, "acceptable_iter", 0)
    problem.x = initial_x
    return problem


end
