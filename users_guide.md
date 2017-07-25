---
layout: default
title: Users' guide
---

[Home](https://chenyx04.github.io/Spectral-Stitching/)  
[Description](https://chenyx04.github.io/Spectral-Stitching/Description)  
[Quick Start](https://chenyx04.github.io/Spectral-Stitching/users_guide)


## Dependencies
The main program uses [NumPy](http://www.numpy.org/) and [SciPy](https://www.scipy.org/) packages, and supports two types of input / output data formats. The evaluation code requires [CVX](http://cvxr.com/cvx/) in MATLAB as well as foreach & doParallel package in R. 


## Usage of Python Main Program


The repo contains a preliminary Python code for the Spectral-Stitching algorithm. It takes as input either the refhap reads file or the contact map file. The output can be either in blocks or a continuous SNP sequence.

### Demo

Two examplanary datasets are provided. Please run the following script to get help
```
python spectral_stitching -h
```
The demo script for running on refhap reads file is 
```
python spectral-stitching.py -r refhapreads_demo.txt  -o out.phased
```

For running on contact map file, use
```
python spectral-stitching.py -c contactmap_demo.csv  -o out.phased
```

If you want a whole continuous snp output instead of blocked one, try parameter -nb.

If no quality score is included or you don't want the quality score to be taken into consideration, use parameter -ns.


### Input File Format 


Each sequenced read is mapped to the reference genomic sequence to obtain the alleles it has at each of the heterozygous sites. Here we represent the major allele by 0 and the minor allele by 1.

The refhap reads format is the same as Refhap algorithm in [Duitama et al. (2010)](http://dl.acm.org/citation.cfm?id=1854802) and Probhap in [V. Kuleshov (2014)](https://www.ncbi.nlm.nih.gov/pubmed/25161223). A demo file is provided in `refhapreads_demo.txt`. The format is as follows.

\# of reads (space) \# of SNPs

\# of SNP contigs covered by read 1 (space) read name (space) start point of first SNP contig (space) alleles in the first SNP contig (space) start point of second SNP contig (space) alleles in the second SNP contig ...

\# of SNP contigs covered by read 2 (space) read name (space) start point of first SNP contig (space) alleles in the first SNP contig (space) start point of second SNP contig (space) alleles in the second SNP contig ...

...



The format of contact map file should be:

\# of rows (space) \# of columns (space) Number of linkages

SNP pos x1 (space) SNP Index y1 (space) \# of reads that indicate x1 and y1 are both major alleles or both minor alleles (If one read indicates x1 and y1 are all 0 or all 1, then this adds 1, else this minuses 1. If the value is less than 0, it indicates that one of x1 and y1 is major allele, and another is minor allele.)

SNP pos x2 (space) SNP Index y2 (space) \# of reads that indicate x2 and y2 are both major alleles or both minor alleles

...

For example, if you have 4 reads, the first read indicates SNP1 and SNP2 are in the same community. The second read indicates SNP1 and SNP2 are in different community. The third read indicate SNP1 and SNP3 are in the same community. The fourth read indicates SNP1 and SNP2 are in different community. Then the contact map should be:

3 3 3

1 1 -1     (The first read +1, the second read -1, the fourth read -1.)

1 3 1      (The third read +1.)

A demo contactmap file is provided in `contactmap_demo.txt`.



### Output File Format

We support two formats of output. By default, the output is in blocks. Each block is printed separately. If you want a contiguous SNP sequence from start to the end ignoring the blocks, please use the -nb or --noblock options.

### Algorithm procedure

When reading in a refhap reads file, we will get the contactmap from all the reads and the quality scores. The blocks are then determined by finding all the connected components from the contact map. The Spectral-Stitching algorithm runs on each block separately and sequentially. In the future release, we will add parallel computing options to our program.


Note: In the new implementation of the Python code, we use the same definition of block as refhap and probhap etc, which is not identical to 10X definition, the one we use for evaluation package. For this version we are showing fewer blocks. The comparison of two different versions is shown in `evaluation\new_implemention_comparison.png`.


## Usage of evaluation package

In the directory `evaluation\`, we provide a demo for running on Chromosome 20 of NA12878 WGS data.

First run `subsample&preprocess.py` for subsampling the coverage and generate adjacent matrix file required for main program.

Then run `pipeline.m` to run Spectral-Stitching on the generated adjacency matrix. 

Finally, run `measure_performance_across_subsamples.R` to get the metrics switch error rate and unphased SNPs. Run `N50.py` to get N50.

If you have an adjacency matrix (or contact map) in the same format as `contactmap_demo.csv`, you can run in MATLAB

```
phased_seq = Spectral_stitching(`contactmap_demo.csv`);
```

to get the phased sequence.

For all the ground truth file for 10X NA12878 WGS data, please go to https://drive.google.com/open?id=0B583FJ-HAtQoY2JFczJhOHlVUnM
