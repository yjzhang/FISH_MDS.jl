#FISH-MDS: Multidimensional scaling (MDS) For Chromosome Structure Inference With Constraints From FISH
## Tutorial By Steven Criscione and Jack Zhang
Email: [steven_criscione@brown.edu](mailto:steven_criscione@brown.edu)
Email: [yue_zhang@alumni.brown.edu](mailto:steven_criscione@alumni.brown.edu)

## Dependencies
FISH_MDS.jl is a Julia package that models 3D structure from Hi-C contact matrices. To use the package please first download Julia and the following libraries.  The companion viewer for the package is Hi-Brow.

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


##FISH_MDS.jl installation:

First, in the command-line open interactive Julia, run 

    julia
    julia> Pkg.clone("https://github.com/yjzhang/FISH_MDS.jl.git")

Then, create a "main.jl" file with the following commands:

    julia> using FISH_MDS
    julia> mds_main()

##Input files

1) The first required file the running main.jl is the normalized Hi-C interaction contact matrix.
In the following examples the file is:  
`HiC_matrix_chr4_condition1.csv`

This file is a comma seperate csv file with no headers or row names.  The file must be normalized contacts and not raw counts.  Raw signal from Hi-C experiments has inherent biases that are disruptive to MDS procedure.
Example lines from HiC_matrix_chr4_condition1.csv:

    104.6821,150.2740,146.4277,0.0000,0.0000,0.0000,...
    995.8,127.0956,137.6080,136.1468,60.2406,53.2425,...
    ...

2) A second useful file to maintain is:  
    `HiC_chr4_bins.bed.` 
    
This is a bed format tab-delimited text file that contains the positions of the bins used by the Hi-C contact matrix.  For example, if the bin size is 200000 (0.2 Mb) then this file will contain:

    chr1    0       200000
    chr1    200000  400000
    chr1    400000  600000

And the file will correspond to the rows of the HiC_matrix_chr4_condition1.csv signal.  Therefore, the number of rows in the Hi-C matrix and the bins file should be identical.  This bed file is not needed by main.jl, but is useful to build optional file 3 and additional tracks for viewing together with the 3D model.


3) The third type input file is optional. This is the file containing measured distances between DNA FISH probes from a 3D DNA FISH experiment typically measured in microns.  In the second running example for main.jl the file is:  
    `chr4_condition1_FISH_probes.txt`

This file looks like:

    3
    7900000 18700000 3.1061
    18700000 26700000 2.8674
    7900000 26700000 1.5315

A space delimited text file. Line-1 is the number of subsequent lines, in this case 3.  Field 1 of line-2 is the bin in the Hi-C matrix for probe-1.  Field 2 of line-1 is the bin in the Hi-C matrix for probe-2.  Field 3 is the distance in microns between probe-1 and probe-2. Line-3 is the interaction between probe-2 and probe-3. Line-4 is the interaction between probe-1 and probe-3.

To find the bin for the probes 1-3, you can intersect the genomic coordinates of the probes (in bed format) with the bins `HiC_chr4_bins.bed.`  (in bed format) using [bedtools](http://bedtools.readthedocs.org/en/latest/#) and the command `bedtools intersect`. In the case a single probe, like probe-1 for instance, overlaps with multiple bins the median bin start position can be used as the position of the probe.

##Examples

The main MDS algorithm used by FISH_MDS.jl is MDS2 from Varoquaux et al. [1].
To run the defualt MDS without constraints run:

`julia main.jl -o chr4_condition1.txt HiC_matrix_chr4_condition1.csv > chr4_condition1.stdout `

To run MDS with provided FISH constraints do:

`julia main.jl -f chr4_sen_probes.txt -o chr4_sen_arm1_test_mds.txt --auto-scale data_fish/combined_Sen_merged_chr4_arm1.csv`

This runs MDS with provided FISH constraints.

`julia main.jl -f data_fish/chr4_qui_probes.txt -o chr4_qui_arm1_test_mds.txt --auto-scale data_fish/combined_Young_Qui_merged_chr4_arm1.csv`


## References

1) Varoquaux, Nelle, Ferhat Ay, William Stafford Noble, and Jean-Philippe Vert. “A Statistical Approach for Inferring the 3D Structure of the Genome.” Bioinformatics 30, no. 12 (June 15, 2014): i26–33. doi:10.1093/bioinformatics/btu268.

2) Lesne, Annick, Julien Riposo, Paul Roger, Axel Cournac, and Julien Mozziconacci. “3D Genome Reconstruction from Chromosomal Contacts.” Nature Methods advance online publication (September 21, 2014). doi:10.1038/nmeth.3104.

