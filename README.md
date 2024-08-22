# geneticmaps

## Table of contents

- [Dependencies](#dependencies)
- [Usage](#usage)
- [Outputs](#outputs)



## Dependencies

### R libraries

- [rqtl](https://rqtl.org/)
- [argparse](https://cran.r-project.org/web/packages/argparse/index.html)

#### Installation using conda

Create a conda environment. 

```
conda create -n r-env
conda activate r-env
```

Then install the packages using conda.

```
conda install -c conda-forge r-essentials r-argparse r-qtl
conda install -c bioconda bcftools (optional)
```

## Usage

### Step1: import.R

```
Rscript import.R -h 
usage: import.R [-h] -i INPUT [-c <CHR> [<CHR> ...]] [-ce <CHR> [<CHR> ...]]
                -pA PARENTA -pB PARENTB --crosstype {bc,f2,riself,risib}

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input VCF file. gzipped file is supported.
  -c <CHR> [<CHR> ...], --chromInclude <CHR> [<CHR> ...]
                        A list of chromosomes to be INcluded in the analysis,
                        separated by space
  -ce <CHR> [<CHR> ...], --chromExclude <CHR> [<CHR> ...]
                        A list of chromosomes to be EXcluded in the analysis,
                        separated by space
  -pA PARENTA, --parentA PARENTA
                        Name of parent A in the vcf file.
  -pB PARENTB, --parentB PARENTB
                        Name of parent B in the vcf file.
  --crosstype {bc,f2,riself,risib}
                        Cross type for the analysis, choose from [bc, f2,
                        riself, risib]
```

### Step2: preprocess.R

```
Rscript preprocess.R -h
usage: preprocess.R [-h] [--filterMissingMarkers]
                    [--filterMissingMarkersThres <FLOAT>] [--filterDupMarkers]
                    [--filterCloseMarkers] [--filterCloseMarkersThres <FLOAT>]
                    [--filterSegregDistMarkers] [--filterMatchingIndividuals]
                    [--filterMatchingIndividualsThres <FLOAT>]

options:
  -h, --help            show this help message and exit
  --filterMissingMarkers
                        Filter out markers with missing data
  --filterMissingMarkersThres <FLOAT>
                        Filter out markers with missing data. For example, if
                        <Thres>=0.05, then any marker with more than 5%
                        missing observations will be removed. (default: 0.05)
  --filterDupMarkers    Filter out markers that have the same genotypes.
  --filterCloseMarkers  Filter out markers that are closer than a specified
                        threshold
  --filterCloseMarkersThres <FLOAT>
                        Identify the largest subset of markers for which no
                        two adjacent markers are separated by less than the
                        specified distance (in cM). (default: 0.5)
  --filterSegregDistMarkers
                        Remove markers that do not pass the segregation test
  --filterMatchingIndividuals
                        Omit individuals with a high proportion of matching
                        genotypes
  --filterMatchingIndividualsThres <FLOAT>
                        Threshold for removing one of a pair of individuals
                        with more than <Thres> proportion of matching
                        genotypes across all markers. (default: 0.9)
```

### Step3: group.R

```
Rscript group.R -h
usage: group.R [-h] [--by {obs,infer}] [--max_rf MAX_RF] [--min_lod MIN_LOD]

options:
  -h, --help         show this help message and exit
  --by {obs,infer}   Specify 'obs' to use predefined groups (observation) or
                     'infer' to use max_rf and min_lod for forming linkage
                     groups (inference)
  --max_rf MAX_RF    Maximum recombination fraction for forming linkage groups
                     (ignored if --group is 'byRef'). (default: 0.25)
  --min_lod MIN_LOD  Minimum LOD score for forming linkage groups (ignored if
                     --group is 'byRef'). (default: 3)

```

P.S. *LOD score*: A statistical estimate of whether two genetic loci are physically near enough to each other (or "linked") on a particular chromosome that they are likely to be inherited together. **A LOD score of 3 or higher is generally understood to mean that two genes are located close to each other on the chromosome.** In terms of significance, a LOD score of 3 means the odds are 1,000:1 that the two genes are linked and therefore inherited together. Also called logarithm of the odds score.

### Step4: order.R

```
Rscript order.R -h      
usage: order.R [-h] [--by {obs,infer}] [-e ERROR_PROB]
               [-f {haldane,kosambi,c-f,morgan}] [-w WINDOW]

options:
  -h, --help            show this help message and exit
  --by {obs,infer}      Specify 'obs' to use predefined orders (observation)
                        but re-estimate the genetic distances or 'infer' to
                        reorder markers (inference)
  -e ERROR_PROB, --error_prob ERROR_PROB
                        Assumed genotyping error rate used in the final
                        estimated map. (default: 0.0001)
  -f {haldane,kosambi,c-f,morgan}, --map_function {haldane,kosambi,c-f,morgan}
                        map function to use. default is 'haldane'.
  -w WINDOW, --window WINDOW
                        window size used to ripple. (default: 3)
```

### Step5: ripple.R

```
Rscript ripple.R -h
usage: ripple.R [-h] [-e ERROR_PROB] [-f {haldane,kosambi,c-f,morgan}]
                [-w WINDOW]

options:
  -h, --help            show this help message and exit
  -e ERROR_PROB, --error_prob ERROR_PROB
                        Assumed genotyping error rate used in the final
                        estimated map. (default 0.0001)
  -f {haldane,kosambi,c-f,morgan}, --map_function {haldane,kosambi,c-f,morgan}
                        map function to use. default is 'haldane'.
  -w WINDOW, --window WINDOW
                        window size used to ripple. (default 3)
```

### Step6: output.R

```
Rscript output.R -h
usage: output.R [-h] [-m]

options:
  -h, --help            show this help message and exit
  -m, --show_marker_names
                        If specified, marker names are included.
```

### Example:

before running these steps, it is highly recommended to use `bcftools` to inspect the raw vcf file.

two commands I find useful are:

```
bcftools query -l [vcf] # to check all individuals' names
```

and

```
bcftools query -f'%CHROM\n' [vcf] | uniq -c # check chromosomes' names
```

#### step1:

then to run this step with parent A named `NM9` and B named `Y158` and crosstype `riself` and without chromosome `Un` :

```
Rscript import.R -i hq.vcf.gz -ce Un -pA NM9 -pB Y158 --crosstype riself
```

```
Copying original file into output/import...
Copy complete.
```

#### step2:

here we only removed duplicate markers and markers that show segregation distortion to simplify.

```
Rscript preprocess.R --filterDupMarkers --filterSegregDistMarkers --filterMatchingIndividuals
```

```
 --Read the following data:
	 299  individuals
	 65364  markers
	 1  phenotypes
 --Cross type: riself 
Warning message:
In convert2riself(cross) : Omitting 529989 genotypes with code==2.
Filtering markers that are duplicates...
Original marker number: 65364, Filtered: 45938, Remaining: 19426
Filtering markers that have segregation distortion...
Original marker number: 19426, Filtered: 1961, Remaining: 17465
Filtering individuals that have more than 90% matching genotypes across all markers...
Original individual number: 299, Filtered: 46, Remaining: 253
```

#### step3: 

```
Rscript group.R --by obs
```

#### step4: 

we recommend re-order the SNPs despite they are already in sequence by reference genome. The rationale is:

Reference genetic maps are often based on a rather small number of individuals. (For example, the original MIT mouse genetic map was based on an intercross with just 46 individuals.) One’s own data often contains many more individuals, and so may produce a more accurate map. The only caveat is that reference genetic maps generally contain a much more dense set of markers, which provides greater ability to detect genotyping errors. Thus reference genetic maps may be based on cleaner genotype data.

```
Rscript order.R --by infer
```

#### step5: (this step can be skipped)

```
Rscript ripple.R -w 2
```

#### step6:

```
Rscript output.R
```

## Outputs

The output files are structured inside the directory `output`.

```
output
├── group
│   ├── map.png
│   ├── mapthis.RDS
│   └── sum.csv
├── import
│   ├── map.png
│   ├── mapthis.RDS
│   ├── org.temp.vcf
│   ├── org.vcf.gz
│   ├── org_map.png
│   └── rqtl.csv
├── order
│   ├── map.png
│   ├── mapthis.RDS
│   └── sum.csv
├── output
│   ├── 1A.png
│   ├── 1B.png
│   ├── 1D.png
│   ...
│   ├── map.png
│   └── sum.csv
└── preprocess
    ├── map.png
    └── mapthis.RDS
```

Each directory contains results of each step.

- `sum.csv`: the chromosomes and position information of SNPs of each step. **The final map information is in** `output/sum.csv` .
- `compare_sum.csv` (if `--by infer` in running `group.R`): the output is a data frame with rows corresponding to the markers and with **two columns**: the initial chromosome assignment and the inferred linkage group. Linkage groups are ordered by the number of markers they contain (from largest to smallest).

+ `map.png`: constructed linkage map plot for every step so far. **The final map plot is in** `output/map.png` .

  ![map](https://github.com/Brycealong/geneticmaps/blob/main/output/output/map.png)



## Citation

- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8496254/
- https://rqtl.org/
- https://rqtl.org/tutorials/geneticmaps.pdf
- https://bioinformaticsworkbook.org/dataAnalysis/GenomeAssembly/GeneticMaps/creating-genetic-maps.html#gsc.tab=0