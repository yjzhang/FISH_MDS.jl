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

    usage: main.jl [-f FISH] [-r RADIUS] [-o OUTPUT] [-i INIT] [--interp]
        [--auto-scale] [--shortest-paths] [-h] filename
    
    positional arguments:
       filename             CSV file containing HiC counts

    optional arguments:
        -f, --fish FISH      FISH constraints file (default: "")
        -r, --radius RADIUS  Radius of chromosome (determined via imaging) (type: Number, default: 1)
        -o, --output OUTPUT  Output file name (default: "")
        -i, --init           Initial starting coords file (in same format as output) (default: "")
        --interp             Flag: use interpolation
        --auto-scale         Flag: automatically infer the scaling factor (requires FISH constraints)
        --shortest-paths     Flag: use shortest paths reconstruction
        -h, --help           show this help message and exit


##FISH_MDS.jl installation:

First, in the command-line, open interactive Julia and run 

    julia
    julia> Pkg.clone("https://github.com/yjzhang/FISH_MDS.jl.git")

Then, create a "main.jl" file with the following commands:

     using FISH_MDS
     mds_main()

You can run this with `julia main.jl [args]`.

Alternatively, in the command line Julia interpreter, run `run_mds(filename, args)`.

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

    chr4    0       200000
    chr4    200000  400000
    chr4    400000  600000

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
1) To run the defualt MDS without constraints run:

`julia main.jl -o chr4_condition1.txt HiC_matrix_chr4_condition1.csv > chr4_condition1.stdout `

The results will be output to `chr4_condition1.txt`. In the results file `chr4_condition1.txt`, the first line states the # of subsequent lines, and the remaining lines contain x,y,z coordinates in space seperated format. Removed indices which are relevant for visualization can be obtained from `chr4_condition1.stdout`, which contains all the details of the MDS optimization run.  

2) To run MDS with provided FISH constraints do:

`julia main.jl -f chr4_condition1_FISH_probes.txt -o chr4_condition1.txt HiC_matrix_chr4_condition1.csv > chr4_condition1.stdout`

This version will fix the relative distances between 3D DNA FISH in the resultant 3D chromosome structure.


3) To run MDS with provided FISH constraints and re-scale to 3D DNA FISH distance run:

`julia main.jl --auto-scale -f chr4_condition1_FISH_probes.txt -o chr4_condition1.txt HiC_matrix_chr4_condition1.csv > chr4_condition1.stdout`

This version will fix the relative distances between 3D DNA FISH in the resultant 3D chromosome structure. It will also scale the resultant 3D structure such that the numerical values of the distances between probes are equal to the distances in microns from the 3D DNA FISH.

4) To run the defualt MDS and re-scale by 3D chromosome paint radius:

`julia main.jl -r 0.7179 -o chr4_condition1.txt HiC_matrix_chr4_condition1.csv > chr4_condition1.stdout`

This version will run default MDS. The resultant 3D structure is re-scaled by the radius calculated from the volume of a chromosome in um^3 measured in a 3D chromosome painting experiment.

5) To run MDS, apply FISH constraints, and re-scale by 3D chromosome paint radius:

`julia main.jl -r 0.7179 -f chr4_condition1_FISH_probes.txt -o chr4_condition1.txt HiC_matrix_chr4_condition1.csv > chr4_condition1.stdout`

This version will fix the relative distances between 3D DNA FISH in the resultant 3D chromosome structure. It will also scale the resultant 3D structure such that the radius is equal to that calculated from the volume in um^3 of a chromosome.  The volume in um^3 is measured from 3D chromosome painting experiments.

6-10) To run a shortest paths computation of wish distances before defualt MDS do:

`julia main.jl --shortest-paths -o chr4_condition1.txt HiC_matrix_chr4_condition1.csv > chr4_condition1.stdout `

All the previous MDS versions can be re-run using precomputation of the wish distances in a shortest paths algorithm as proposed by Lesne et. al. [2]. 

11) To run an MDS computation using an initialized set of wish distances run:

`julia main.jl -i chr4_control1.txt -o chr4_condition1.txt HiC_matrix_chr4_condition1.csv > chr4_condition1.stdout `

This version will try to run MDS computation from an initial set of wish distances provided by a control rather than hypothesizing wish distances.  This version of the command prevents artefacts such as global reflections of the solution for a condition with respect to a control sample. 

##Removing outlier coordinates

We provide a utility script to remove outlier coordinates from the resultant structural solution.  This script identifies outliers from the distribution of distances between coordinates and replaces them with the midpoint of adjacent coordinates.  
`Rscript remove_outliers.R chr4_condition1.txt `

This script creates the output file `fix.chr4_condition1.txt ` with the outlier coordinates removed.  The utility script is located within the `src` folder which can be added to the `$PATH` variable to be used from anywhere. 

##Re-binning data for visualization

We also provide a second utility script for re-binning genomic signal data (such as ChIP-seq) for compatability with Hi-C binning.  This is recommended for only broad signal profiles, such as certain histone modifications. The script requires [bedtools](http://bedtools.readthedocs.org/en/latest/#), the python module [pybedtools] (https://pythonhosted.org/pybedtools/), and for bigWig format requires the additional binary utility [bigWigToBedGraph] (http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/bigWigToBedGraph) in the `$PATH` variable. To add a script to  `$PATH`  you can do the following (make sure script has executable priveledges):
```
working_path=/path/to/script_folder
export PATH=${working_path}:$PATH

```

To run this utility script run the command:
`python rebin_bedgraph.py input_signal bins resolution outname`

an example would be:
`python rebin_bedgraph.py H3K27me3_signal.bedGraph HiC_chr4_bins.bed 200000 H3K27me3_rebinned.bedGraph`

For more information see detailed help page:
`python rebin_bedgraph.py --help`

A bedgraph file is a tab seperated genomic signal file without a header:
```chr     start   end     signal```
for example:
```
chr4   50000   75000   1.56248
chr4   75000   100000  2.21352
chr4   100000  125000  1.7006
```
The tool will re-bin the signal, and return the re-binned data in the following format with a header:
```chr   start   end   signal   log2.signal   rescaled.signal   trimmed.signal   log2trimmed.signal```
Where ```signal``` is the mean signal for Hi-C bins, ```log2.signal``` is the log2 transform of that mean signal.  The ```rescaled.signal``` is normalized mean signal to range (0-1) via normalized = (x-min(x))/(max(x)-min(x)).  The ```trimmed.signal``` is the [Winsorisation] (https://en.wikipedia.org/wiki/Winsorising) transform of the data removing the bottom and top 5th percentile for the mean signal.  The ```log2trimmed.signal``` is the [Winsorisation] (https://en.wikipedia.org/wiki/Winsorising) transform of the log2 data removing the 5th percentile for the mean signal.  

## References

1) Varoquaux, Nelle, Ferhat Ay, William Stafford Noble, and Jean-Philippe Vert. “A Statistical Approach for Inferring the 3D Structure of the Genome.” Bioinformatics 30, no. 12 (June 15, 2014): i26–33. doi:10.1093/bioinformatics/btu268.

2) Lesne, Annick, Julien Riposo, Paul Roger, Axel Cournac, and Julien Mozziconacci. “3D Genome Reconstruction from Chromosomal Contacts.” Nature Methods advance online publication (September 21, 2014). doi:10.1038/nmeth.3104.

