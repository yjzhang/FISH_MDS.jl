include("mds_approximate_fish.jl")
include("output.jl")
include("interpolate.jl")

filename = "data_fish/combined_Sen_merged_chr4_arm1.csv"
data = load_data(filename, scale=1, exponent=1/3)
#mds = MDSProblem(data, 1.0)#, lads, nads, 0.1)
mds = new_mds(data, 1.0, "data_fish/chr4_sen_probes.txt")
mds = remove_infs(mds)
println(filename)
#println(mds)
ipopt_problem = make_ipopt_problem(mds, radius_constraint=true,
    fish_constraint=true)
solveProblem(ipopt_problem)
coords = ipopt_problem.x
coords = reshape(coords, 3, int(length(coords)/3))'
interp_coords = interp_3d(coords, 10)
output_txt(coords, "chr4_sen_arm1_with_constraint.txt")
output_txt(interp_coords, "interp_chr4_sen_arm1_with_constraint.txt")

ipopt_problem = make_ipopt_problem(mds, radius_constraint=true,
    fish_constraint=false)
solveProblem(ipopt_problem)
coords = ipopt_problem.x
coords = reshape(coords, 3, int(length(coords)/3))'
output_txt(coords, "chr4_sen_arm1_no_constraint.txt")

filename = "data_fish/combined_Young_Qui_merged_chr4_arm1.csv"
data = load_data(filename, scale=1, exponent=1/3)
#mds = MDSProblem(data, 1.0)#, lads, nads, 0.1)
mds = new_mds(data, 1.0, "data_fish/chr4_qui_probes.txt")
mds = remove_infs(mds)
println(filename)
#println(mds)
ipopt_problem = make_ipopt_problem(mds, radius_constraint=true,
    fish_constraint=true)
solveProblem(ipopt_problem)
coords = ipopt_problem.x
coords = reshape(coords, 3, int(length(coords)/3))'
interp_coords = interp_3d(coords, 10)
output_txt(coords, "chr4_qui_arm1_with_constraint.txt")
output_txt(interp_coords, "interp_chr4_qui_arm1_with_constraint.txt")


filename = "data_fish/combined_You_Pro_merged_chr4_arm1.csv"
data = load_data(filename, scale=1, exponent=1/3)
mds = MDSProblem(data, 1.0)#, lads, nads, 0.1)
#mds = new_mds(data, 1.0, "data_fish/chr4_qui_probes.txt")
mds = remove_infs(mds)
println(filename)
#println(mds)
ipopt_problem = make_ipopt_problem(mds, radius_constraint=true,
    fish_constraint=false)
solveProblem(ipopt_problem)
coords = ipopt_problem.x
coords = reshape(coords, 3, int(length(coords)/3))'
interp_coords = interp_3d(coords, 10)
output_txt(coords, "chr4_pro_arm1_no_constraint.txt")
output_txt(interp_coords, "interp_chr4_pro_arm1_no_constraint.txt")


filename = "../fish/1M/combined_Sen_merged_chr4_arm1.csv"
data = load_data(filename, scale=1, exponent=1/3)
mds = MDSProblem(data, 1.0)#, lads, nads, 0.1)
#mds = new_mds(data, 1.0, "data_fish/chr4_sen_probes.txt")
mds = remove_infs(mds)
println(filename)
#println(mds)
ipopt_problem = make_ipopt_problem(mds, radius_constraint=true,
    fish_constraint=false)
solveProblem(ipopt_problem)
coords = ipopt_problem.x
coords = reshape(coords, 3, int(length(coords)/3))'
interp_coords = interp_3d(coords, 10)
output_txt(coords, "1M_chr4_sen_arm1_with_constraint.txt")
output_txt(interp_coords, "1M_interp_chr4_sen_arm1_with_constraint.txt")
