# basically, writing stuff out in a nice format so that other things can read them

SCALE_FACTOR = 100

function output_pdb(coords::Array{Float64, 2}, out_file::AbstractString)
    # coords: each row is a coord
    n, dim = size(coords)
    out = open(out_file, "w")
    for i in 1:n
        c = coords[i, :]
        write(out, @sprintf("ATOM  "))
        write(out, @sprintf("%5d", i))
        write(out, @sprintf(" "))
        write(out, @sprintf("N   ")) # atom name
        write(out, @sprintf(" "))
        write(out, @sprintf("NOD"))
        write(out, @sprintf(" "))
        write(out, @sprintf("%c", 'A'))
        write(out, @sprintf("    "))
        write(out, @sprintf(" "))
        write(out, @sprintf("   "))
        write(out, @sprintf("%8.3f%8.3f%8.3f", (c[1]+1)*SCALE_FACTOR, (c[2]+1)*SCALE_FACTOR, (c[3]+1)*SCALE_FACTOR))
        write(out, @sprintf("%6.2f", 1.0))
        write(out, @sprintf("%6.2f", 50.0))
        write(out, @sprintf(" "))
        write(out, @sprintf("  "))
        write(out, @sprintf("\n"))
    end
    close(out)
end

function output_txt(coords::Array{Float64, 2}, out_file::AbstractString)
    # coords: each row is a coord
    # outputs a space-separated text file
    n, dim = size(coords)
    out = open(out_file, "w")
    write(out, @sprintf("%d\n", n))
    for i in 1:n
        c = coords[i,:]
        write(out, @sprintf("%8.6f ", c[1]))
        write(out, @sprintf("%8.6f ", c[2]))
        write(out, @sprintf("%8.6f ", c[3]))
        write(out, "\n")
    end
    close(out)
end
