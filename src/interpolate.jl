using Grid

function interp_3d(data::Array{Float64, 2}, num_points::Int)
    # interpolates with num_points per original point
    interp_type = InterpQuadratic
    y1 = InterpGrid(data[:,1], BCnil, interp_type)
    y2 = InterpGrid(data[:,2], BCnil, interp_type)
    y3 = InterpGrid(data[:,3], BCnil, interp_type)
    output = zeros((size(data, 1)-1)*num_points + 1, 3)
    for i in 1:(size(data, 1)-1)
        index = 1
        for j in range(0, 1/num_points, num_points)
            output[(i-1)*num_points + index, :] = [y1[i+j], y2[i+j], y3[i+j]]
            index += 1
        end
    end
    output[end, :] = data[end, :]
    return output
end

function centroid(data::Array{Float64, 2})
    # Returns the centroid of the points - a 3x1 vector.
    return mean(data, 1)
end

function center_data(data::Array{Float64, 2})
    # centers the data such that the centroid is at 0.
    centroid = mean(data, 1)
    return data .- centroid
end
