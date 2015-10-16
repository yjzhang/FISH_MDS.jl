
module FISH_MDS

    using Ipopt
    using ArgParse
    using Grid

    # exports...
    
    export run_mds
    export mds_main

    # things to run w/o command line?
    export load_data
    export new_mds
    export make_ipopt_problem
    export solveProblem
    export remove_infs
    export get_removed_indices
    export output_txt
    export interp_3d

    include("graph.jl")
    include("mds_metric.jl")
    include("mds_approximate_fish.jl")
    include("output.jl")
    include("interpolate.jl")
    include("mds_main.jl")


end # module
