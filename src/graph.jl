type GraphNode
    index::Int
end

type Graph
    nodes::Dict{Int, GraphNode}
    edges::Dict{Int, Array{Int, 1}}
    lengths::Dict{Tuple{Int, Int}, Real}
    max_index::Int
end

function init_graph()
    return Graph(Dict{Int, GraphNode}(), Dict{Int, Array{Int, 1}}(),
        Dict{Tuple{Int, Int}, Real}(), 1)
end

function num_nodes(g::Graph)
    return length(g.nodes)
end

function add_node!(g::Graph, index::Int)
    # Adds a node given its neighbors and edge lengths.
    new_node = GraphNode(index)
    g.nodes[index] = new_node
    g.edges[index] = Array(Int, 0)
    return index 
end

function add_node!(g::Graph, index::Int, n::Array{Int, 1}, l::Array)
    # Adds a node given its neighbors and edge lengths.
    new_node = GraphNode(index)
    if index in keys(g.nodes)
        println("Warning: index already exists")
    end
    g.nodes[index] = new_node
    g.max_index += 1
    g.edges[index] = n
    for (node, length) in zip(n, l)
        push!(g.edges[node], index)
        g.lengths[(node, index)] = length
        g.lengths[(index, node)] = length
    end
    return index 
end

function add_edge!(g::Graph, i1::Int, i2::Int, d::Real)
    push!(g.edges[i1], i2)
    push!(g.edges[i2], i1)
    g.lengths[(i1, i2)] = d
    g.lengths[(i2, i1)] = d
end

function all_shortest_paths(g::Graph)
    # returns the all-source shortest paths as a NxN matrix
    # Using Floyd-Warshall
    n = num_nodes(g)
    dist = fill(Inf, n, n)
    for (i1, i2) in keys(g.lengths)
        dist[i1, i2] = g.lengths[i1, i2]
    end
    for i in keys(g.nodes)
        dist[i, i] = 0
    end
    for k in range(1, n)
        for i in range(1, n)
            for j in range(1, n)
                if dist[i,j] > dist[i,k] + dist[k,j]
                    dist[i,j] = dist[i,k] + dist[k,j]
                end
            end
        end
    end
    return dist
end
