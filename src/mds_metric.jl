include("graph.jl")

function contact_freq_to_distance_matrix(contact_file::AbstractString; remove_inf=true, scale=1, exponent=1/3)
    # turns a csv file of contact frequences into a distance matrix
    contact_data = readcsv(contact_file, Float64)
    # TODO: use the graphical approach... distances are 1/freq, use shortest
    # paths algorithm to get distances that don't exist

    # arbitrarily scale distances by some scale
    distances = scale./(contact_data).^exponent
    dist = distances
    # removing zero points
    if remove_inf
        dist_index = 1
        while dist_index <= size(dist)[1]
            if mapreduce(x->(isinf(x) || x==0), &, dist[dist_index, 1:end-1])
                dist = [dist[1:dist_index-1, :], dist[dist_index+1:end, :]]
                dist = [dist[:, 1:dist_index-1] dist[:,dist_index+1:end]]
            else
                dist_index += 1
            end
        end
    else
         
    end
    distances = dist 
    n, n2 = size(distances)
    g = init_graph()
    for i in range(1, n)
        add_node!(g, i)
    end
    for i in range(1, n-1)
        for i2 in (i+1):n
            if distances[i, i2] != Inf
                add_edge!(g, i, i2, distances[i, i2])
            else
                #distances[i, i2] = 1.0
                #add_edge!(g, i, i2, distances[i, i2])
            end
        end
    end
    # done with constructing graph, now to find the distances
    dist = all_shortest_paths(g)
    # TODO: what if there's an "Inf" in the distance? replace by 1?
    # this is an incredibly hacky way of doing things, and it might mess
    # everything up entirely.
    # ... and it does seem to mess up everything entirely.
    # how else would this be done?
    #dist[dist.==Inf] = 1
    return dist
end

function mds_metric(dist::Array{Float64, 2}; dims = 3)
    # 'dist' is a matrix of distances (should be symmetric)
    # metric MDS - assumes that the data is metric
    # TODO: do some crazy matrix operations
    N, N = size(dist)
    d = zeros(N)
    b = 1/N^2*sum([sum([dist[j,k]^2 for k in j:N]) for j in range(1,N)])
    for i in range(1, N)
        d[i] = 1/N * sum([dist[i, j]^2 for j in range(1, N)]) - b
    end
    M = zeros(N, N)
    for i in range(1, N)
        for j in range(1, N)
            M[i,j] = 1/2*(d[i] + d[j] - dist[i,j]^2)
        end
    end
    lambda, E = eig(M)
    # largest 3 eigenvalues
    v = E[:, end-dims+1:end] .* sqrt(lambda[end-dims+1:end])' 
    return fliplr(v[:,1:dims]) #[v[:,3] v[:,2] v[:,1]]
end

function test()
    m = "../test_data/combined_You_Pro_merged_chr18.csv"
    dist = contact_freq_to_distance_matrix(m)
    return dist, mds_metric(dist)
end
