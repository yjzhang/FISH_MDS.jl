using FISH_MDS

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
output_txt(coords, "test_out/chr4_sen_arm1_with_constraint.txt")
output_txt(interp_coords, "test_out/interp_chr4_sen_arm1_with_constraint.txt")


# tests: no constraint, with constraint, no interp,
# with interp, etc.

filename = "data_fish/combined_Sen_merged_chr4_arm1.csv"

println("no interp, with constraint")

run_mds(filename, interp=false, 
    output_name = "test_out/sen_chr4_constrained_interp.txt", 
    constraint="data_fish/chr4_sen_probes.txt")

println("no interp, with constraint, with autoscale")

run_mds(filename, interp=false, auto_scale=true, 
    output_name = "test_out/sen_chr4_constrained_interp_autoscale.txt", 
    constraint="data_fish/chr4_sen_probes.txt")

println("interp, with constraint")

run_mds(filename, interp=true, 
    output_name = "test_out/sen_chr4_constrained_interp.txt", 
    constraint="data_fish/chr4_sen_probes.txt")

println("no interp, no constraint")

run_mds(filename, interp=false,
    output_name = "test_out/sen_chr4_default.txt")

println("no interp, no constraint, fixed start")

run_mds(filename, interp=false,
    output_name = "test_out/sen_chr4_default_2.txt", 
    starting_points_file="test_out/sen_chr4_default.txt")
# new file
filename = "data_fish/combined_Young_Qui_merged_chr18.csv"

run_mds(filename, output_name = "test_out/qui_chr18_default.txt")

