include("mds_approximate.jl")
include("mds_metric.jl")
include("output.jl")

function basic_test()
    filename = "test/combined_You_Pro_merged_chr19.csv"
    lads = "test/lads_chr18.txt"
    nads = "test/nads_chr18.txt"
    data = load_data(filename, scale=1, exponent=1)
    mds = MDSProblem(data, 1.0)
    mds = remove_infs(mds)
    ipopt_problem = make_ipopt_problem(mds, radius_constraint=false, 
        lad_constraint=true, nad_constraint=false)
    solveProblem(ipopt_problem)
    x = ipopt_problem.x
    coords = reshape(x, 3, int(length(x)/3))'
    output_txt(coords, "basic_mds_test_19.txt")
end

filename = "../fish/combined_Sen_merged_chr4_arm1.csv"
#filename = "../40kb_test_data/combined_You_Pro_merged_chr18.csv"
#filename = "..test_data/senescent/combined_Sen_merged_chr18.csv"
data = load_data(filename, scale=1, exponent=1)
#data = contact_freq_to_distance_matrix(filename, remove_inf=false, scale=2, exponent=1)
#println(data)
mds = MDSProblem(data, 1.0)#, lads, nads, 0.1)
mds = remove_infs(mds)
#println(mds.dist)
#println(mds.lads)
ipopt_problem = make_ipopt_problem(mds, radius_constraint=true, 
    lad_constraint=false, nad_constraint=false)
solveProblem(ipopt_problem)
coords = ipopt_problem.x
#x = ipopt_problem.x
#println("nuc center: ", x[end-2], " ", x[end-1]," ",  x[end])
coords = reshape(coords, 3, int(length(coords)/3))'
output_txt(coords, "chr4_sen_arm1.txt")
#basic_test()
