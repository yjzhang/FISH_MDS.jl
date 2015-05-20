#FISH-MDS: MDS For Chromosome Structure Inference With Constraints From FISH
## Tutorial By Steven Criscione and Jack Zhang
Email: [steven_criscione@brown.edu](mailto:steven_criscione@brown.edu)
Email: [yue_zhang@alumni.brown.edu](mailto:steven_criscione@alumni.brown.edu)

## Dependencies

[Julia](http://julialang.org/)

Julia libraries:

[Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)  
[Grid.jl](https://github.com/timholy/Grid.jl) (for 3D interpolation)  
[ArgParse.jl](https://github.com/carlobaldassi/ArgParse.jl) (for the command-line interface)  

Visualization: add link to Hi-Brow

## Usage

    usage: main.jl [-f FISH] [-r RADIUS] [-o OUTPUT] [--interp]
        [--auto-scale] [--shortest-paths] [-h] filename
    
    positional arguments:
       filename             CSV file containing HiC counts

    optional arguments:
        -f, --fish FISH      FISH constraints file (default: "")
        -r, --radius RADIUS  Radius of chromosome (determined via imaging) (type: Number, default: 1)
        -o, --output OUTPUT  Output file name (default: "")
        --interp             Flag: use interpolation
        --auto-scale         Flag: automatically infer the scaling factor (requires FISH constraints)
        --shortest-paths     Flag: use shortest paths reconstruction
        -h, --help           show this help message and exit







## Examples

`julia main.jl -o chr4_pro_arm1_test_mds.txt data_fish/combined_You_Pro_merged_chr4_arm1.csv`

This runs MDS without any constraints.

`julia main.jl -f data_fish/chr4_sen_probes.txt -o chr4_sen_arm1_test_mds.txt --auto-scale data_fish/combined_Sen_merged_chr4_arm1.csv`

This runs MDS with provided FISH constraints.

`julia main.jl -f data_fish/chr4_qui_probes.txt -o chr4_qui_arm1_test_mds.txt --auto-scale data_fish/combined_Young_Qui_merged_chr4_arm1.csv`

This runs MDS with provided FISH constraints.

Visualization: `./mds_vis chr4_qui_arm1_test_mds.txt` 

(after running the previous command for qui data)

## References

Varoquaux, Nelle, Ferhat Ay, William Stafford Noble, and Jean-Philippe Vert. “A Statistical Approach for Inferring the 3D Structure of the Genome.” Bioinformatics 30, no. 12 (June 15, 2014): i26–33. doi:10.1093/bioinformatics/btu268.

and

Lesne, Annick, Julien Riposo, Paul Roger, Axel Cournac, and Julien Mozziconacci. “3D Genome Reconstruction from Chromosomal Contacts.” Nature Methods advance online publication (September 21, 2014). doi:10.1038/nmeth.3104.

