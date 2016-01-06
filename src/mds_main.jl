#!/usr/bin/env julia

function get_distance(coords::Array{Float64, 2}, c1::Int, c2::Int)
    # Returns the distance between two points
    return norm(coords[c1, :] - coords[c2, :])
end

function validate_fish_constraints(coords::Array{Float64, 2}, mds::MDSProblem)
    # checks that the FISH constraints are sort of close...
    # TODO: not yet implemented
    c = size(coords, 1)/size(mds.dist, 1)
end

function autoscale(coords::Array{Float64, 2}, mds::MDSProblem)
    # finds the scale that makes the FISH inferred distances correct
    d1 = mds.c1_c2_distance
    return d1/get_distance(coords, mds.control_locus_1, mds.control_locus_2)
end

function radius_scale(coords::Array{Float64, 2}, radius::Real)
    # finds the scale such that the radius of the chromosome becomes the
    # given radius.
    coords = center_data(coords)
    max_radius = maximum([norm(coords[i, :]) for i in 1:size(coords, 1)])
    return radius/max_radius
end

function run_mds(filename::AbstractString; scale::Number=1, constraint::AbstractString="",
            interp::Bool=true, exponent::Number=1/3, output_name::AbstractString="",
            auto_scale::Bool=false, shortest_paths::Bool=false,
            starting_points_file::AbstractString="")
    # loading data
    data = load_data(filename, scale=1, exponent=exponent)
    if shortest_paths
        data = contact_freq_to_distance_matrix(filename, scale=1, 
                exponent=exponent, remove_inf=false)
    end
    mds = MDSProblem(data)
    # auto-creating output file name
    if output_name == ""
        prefix = split(filename, ".")[1]
        output_name = string(prefix, "_mds_coords.txt")
    end
    # initializing MDS problem
    mds = MDSProblem(data)
    zero_indices = get_removed_indices(mds)
    mds = remove_infs(mds)
    println("Running MDS on ", filename, "\n")
    println("Removed indices: ", zero_indices)
    println()
    ipopt_problem = make_ipopt_problem(mds, radius_constraint=false,
            fish_constraint=false)
    # dealing with constraints
    if constraint != ""
        mds = new_mds(data, 1.0, constraint)
        mds = remove_infs(mds)

        print_constraints(mds)
        println()
        ipopt_problem = make_ipopt_problem(mds, radius_constraint=false,
                fish_constraint=true)
    end
    # initial point
    if length(starting_points_file) > 0
        start_points = float64(readdlm(starting_points_file)[2:end,:])
        starting_x = reshape(start_points', length(start_points))
        if length(starting_x) != length(ipopt_problem.x)
            l1 = length(starting_x)
            l2 = length(ipopt_problem.x)
            if l1 < l2
                println("Error: not enough starting coords. Filling with random starting points.")
                # TODO: error?
                starting_x = cat(1, starting_x, ipopt_problem.x[l1+1:end])
            else
                println("Error: too many starting coords. Truncating.")
                println(l1)
                println(l2)
                starting_x = starting_x[1:l2]
            end
        else
            ipopt_problem.x = starting_x
        end
    end
    # solve problem
    solveProblem(ipopt_problem)
    coords = reshape(ipopt_problem.x, 3, round(Int, length(ipopt_problem.x)/3))'
    # provided scale
    if !auto_scale && scale != 1
        scale = radius_scale(coords, scale);
    end
    # autoscale
    if auto_scale && scale == 1 && constraint != ""
        scale = autoscale(coords, mds)
    end
    println("\nInferred scaling factor: ", scale, "\n")
    coords = center_data(coords)
    coords = scale * coords
    # interpolation
    if interp
        coords = interp_3d(coords, 10)
    end
    # file output
    println("\nWriting to file ", output_name)
    output_txt(coords, output_name)
    return coords
end


function parse_cl()
    # uses ArgParse to parse the command line
    s = ArgParseSettings()

    @add_arg_table s begin
        "filename"
            help = "CSV file containing HiC counts"
            required = true
        "--fish", "-f"
            help = "FISH constraints file"
            default = ""
        "--radius", "-r"
            help = "Radius of chromosome (determined via imaging)"
            arg_type = Float64
            default = 1.0
        "--output", "-o"
            help = "Output file name"
            default = ""
        "--init", "-i"
            help = "Starting coords file for optimization"
            default = ""
        "--interp"
            help = "Flag: use interpolation"
            action = :store_true
        "--auto-scale"
            help = "Flag: automatically infer the scaling factor (requires FISH constraints)"
            action = :store_true
        "--shortest-paths"
            help = "Flag: use shortest paths reconstruction"
            action = :store_true
    end

    return parse_args(s)
end

function mds_main()
    args = parse_cl()
    run_mds(args["filename"], 
            constraint = args["fish"], 
            scale = args["radius"], 
            output_name = args["output"],
            auto_scale = args["auto-scale"],
            interp = args["interp"],
            shortest_paths = args["shortest-paths"],
            starting_points_file = args["init"])
end
