
module FISH_MDS

    using Ipopt
    using ArgParse
    using Grid

    # exports...
    
    export run_mds
    export mds_main

    include("graph.jl")
    include("mds_metric.jl")
    include("mds_approximate_fish.jl")
    include("output.jl")
    include("interpolate.jl")
    include("main.jl")


end # module
