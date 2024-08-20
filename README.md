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
```

## Usage

### Step1: import.R

```
Rscript import.R -h 
usage: import.R [-h] -i INPUT [-c <CHR> [<CHR> ...]] [-ce <CHR> [<CHR> ...]]
                -pA PARENTA -pB PARENTB

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
```

### Step2: preprocess.R

```
Rscript preprocess.R -h
usage: preprocess.R [-h] --crosstype {bc,f2,riself,risib}
                    [--filterMissingMarkers]
                    [--filterMissingMarkersThres <FLOAT>] [--filterDupMarkers]
                    [--filterCloseMarkers] [--filterCloseMarkersThres <FLOAT>]
                    [--filterSegregDistMarkers] [--filterMatchingIndividuals]
                    [--filterMatchingIndividualsThres <FLOAT>]

options:
  -h, --help            show this help message and exit
  --crosstype {bc,f2,riself,risib}
                        Cross type for the analysis, choose from [bc, f2,
                        riself, risib]
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

#### step1:

```
$ Rscript import.R -i hq.vcf.gz -ce Un -pA NM9 -pB Y158
Copying original file into output/import...
Copy complete.
```

#### step2:

```
$ Rscript preprocess.R --crosstype riself --filterMissingMarkers --filterDupMarkers --filterCloseMarkers --filterCloseMarkersThres 5 --filterSegregDistMarkers --filterMatchingIndividuals
 --Read the following data:
	 299  individuals
	 65364  markers
	 1  phenotypes
 --Cross type: riself 
Warning message:
In convert2riself(cross) : Omitting 529989 genotypes with code==2.
Filtering markers that have more than 5% missing observations...
Original marker number: 65364, Filtered: 5616, Remaining: 59748
Filtering markers that are duplicates...
Original marker number: 59748, Filtered: 43583, Remaining: 16165
Identifying the largest subset of markers where no two adjacent markers are separated by less than 5 cM...
Original marker number: 16165, Filtered: 14650, Remaining: 1515
Filtering markers that have segregation distortion...
Original marker number: 1515, Filtered: 239, Remaining: 1276
Filtering individuals that have more than 90% matching genotypes across all markers...
Original individual number: 299, Filtered: 46, Remaining: 253
        n.mar  length ave.spacing max.spacing
1A         32   597.3        19.3       286.5
1B         61   694.8        11.6       118.3
1D         39   481.7        12.7       110.4
2A         78   780.3        10.1        54.2
2B        101   806.9         8.1        91.6
2D         79   651.3         8.3        37.7
3A         57   751.4        13.4       398.9
3B         56   843.4        15.3       114.7
3D         47   610.7        13.3       292.0
4A         91   747.9         8.3        45.3
4B         74   669.7         9.2        50.1
4D         25   508.6        21.2       224.0
5A         96   676.7         7.1        27.3
5B         96   711.6         7.5        44.8
5D         54   544.4        10.3       164.4
6A         36   615.6        17.6        58.2
6B         83   708.9         8.6        54.3
6D         39   481.8        12.7        85.3
7A         33   740.8        23.2       546.2
7B         38   759.1        20.5       457.0
7D         61   638.8        10.6        94.4
overall  1276 14021.7        11.2       546.2

```

#### step3:

```
$ Rscript group.R --by obs
null device 
          1 
Using the input groups...
        n.mar  length ave.spacing max.spacing
1A         32   597.3        19.3       286.5
1B         61   694.8        11.6       118.3
1D         39   481.7        12.7       110.4
2A         78   780.3        10.1        54.2
2B        101   806.9         8.1        91.6
2D         79   651.3         8.3        37.7
3A         57   751.4        13.4       398.9
3B         56   843.4        15.3       114.7
3D         47   610.7        13.3       292.0
4A         91   747.9         8.3        45.3
4B         74   669.7         9.2        50.1
4D         25   508.6        21.2       224.0
5A         96   676.7         7.1        27.3
5B         96   711.6         7.5        44.8
5D         54   544.4        10.3       164.4
6A         36   615.6        17.6        58.2
6B         83   708.9         8.6        54.3
6D         39   481.8        12.7        85.3
7A         33   740.8        23.2       546.2
7B         38   759.1        20.5       457.0
7D         61   638.8        10.6        94.4
overall  1276 14021.7        11.2       546.2
```

#### step4:

```
$ Rscript order.R --by obs
Using the input orders...
        n.mar  length ave.spacing max.spacing
1A         32   208.3         6.7        34.1
1B         61  1202.6        20.0       966.8
1D         39   277.6         7.3        55.0
2A         78   322.3         4.2        33.1
2B        101   294.8         2.9        19.1
2D         79   391.1         5.0        33.1
3A         57   286.6         5.1        29.8
3B         56  2198.2        40.0       966.8
3D         47   346.4         7.5        46.6
4A         91   537.7         6.0       135.8
4B         74   198.6         2.7        16.6
4D         25   182.1         7.6        37.2
5A         96   332.6         3.5        22.9
5B         96   342.7         3.6        39.4
5D         54   331.6         6.3        19.8
6A         36  1332.3        38.1       966.8
6B         83   245.5         3.0        25.4
6D         39   251.1         6.6        33.6
7A         33   276.8         8.6        22.5
7B         38   241.6         6.5        21.4
7D         61   414.7         6.9        36.4
overall  1276 10215.5         8.1       966.8
```

#### step5: (this step can be skipped)

```
$ Rscript ripple.R -w 2
Ripple the order using window size 2
   32 total orders
    --Order 5 
    --Order 10 
    --Order 15 
    --Order 20 
    --Order 25 
    --Order 30 
   61 total orders
    --Order 10 
    --Order 20 
    --Order 30 
    --Order 40 
    --Order 50 
    --Order 60 
   39 total orders
    --Order 5 
    --Order 10 
    --Order 15 
    --Order 20 
...
```

#### step6:

```
Rscript output.R
null device 
          1 

```

## Outputs

Files inside `OUT_DIR` are like below.

```
output
├── group
│   ├── compare_sum.csv
│   ├── mapthis.RDS
│   ├── rf-vs-LOD.png
│   └── sum.csv
├── import
│   ├── org.temp.vcf
│   ├── org.vcf.gz
│   └── rqtl.csv
├── order
│   ├── mapthis.RDS
│   └── sum.csv
├── output
│   └── map.png
├── preprocess
│   └── mapthis.RDS
└── ripple
│   ├── mapthis.RDS
└── └── sum.csv
```

Each directory contains results of each step.

+ `SNPindex.tsv` : columns in this order.
  + **CHROM** - The chromosome this SNP is in 
  + **POS** - The position on the chromosome in nt 
  + **REF** - The reference allele at that position 
  + **ALT** - The alternate allele 
  + **DP.HIGH** - The read depth at that position in the high bulk 
  + **AD_REF.HIGH** - The allele depth of the reference allele in the high bulk 
  + **AD_ALT.HIGH** - The alternative allele depth in the high bulk  
  + **SNPindex.HIGH** - The calculated SNP-index for the high bulk 
  + Same as above for the low bulk 
  + **REF_FRQ** - The reference allele frequency as defined above 
  + **deltaSNP** - The $\Delta$(SNP-index) as defined above
+ `SNPindex.filt.tsv` : SNPs filtered with user-specified or default thresholds. One column `tricubeDeltaSNP` is added, which represents the smoothed deltaSNP values.
+ `allchr.png` : delta SNP index for all chromosomes
  - **dots** : variant
  - **<span style="color: red; ">RED line</span>** : smoothed delta SNP-index
    ![allchr](https://github.com/Brycealong/QTL-analysis/blob/main/images/allchr.png)
+ `chr1.png` : delta SNP index for one chromosome. Same for other chromosomes.
  ![6a](https://github.com/Brycealong/QTL-analysis/blob/main/images/6A.png)
+ `distribution.png`: distribution of reference allele frequency, read depths of each sample and SNP index of each sample. Adjust your thresholds using this graph.
  ![dis](https://github.com/Brycealong/QTL-analysis/blob/main/images/distribution.png)
+ `analysis.log`: log how many SNPs are filtered out on each parameter.



## Citation

- Hiroki Takagi, Akira Abe, Kentaro Yoshida, Shunichi Kosugi, Satoshi Natsume, Chikako Mitsuoka, Aiko Uemura, Hiroe Utsushi, Muluneh Tamiru, Shohei Takuno, Hideki Innan, Liliana M. Cano, Sophien Kamoun, Ryohei Terauchi (2013).  [QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations](https://doi.org/10.1111/tpj.12105). Plant journal 74:174-183.
- Mansfeld, B.N. and Grumet, R. (2018), QTLseqr: An R Package for Bulk Segregant Analysis with Next-Generation Sequencing. The Plant Genome, 11: 180006. https://doi.org/10.3835/plantgenome2018.01.0006