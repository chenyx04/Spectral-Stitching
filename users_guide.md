---
layout: default
title: Users' guide
---

[Home](https://chenyx04.github.io/Spectral-Stitching/)

## Usage of Python Main Program


We`ve uploaded a preliminary python version for spectral stitching algorithm. It takes in either the refhap reads file or contact map file. And the output can be either in blocks or a continuous SNP sequence.

### Demo

We provide two examplanary dataset, run the following script to get help
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

If you want a whole continuous snp output instead of blocked one. Try parameter -nb.



### Input File Format 


The refhap reads format is the same as described in the paper by Duitama et al. (2010) and the documentation of RefHap algorithm. A demo file is provided in `refhapreads_demo.txt`.

The format of contact map file should be:

\# of rows (space) \# of columns (space) Number of linkages

SNP pos x1 (space) SNP Index y1 (space) \# of reads that indicate x1 and y1 belong to one community (If one read indicates x1 and y1 are all 0 or all 1, then this adds 1, else this minuses 1)

SNP pos x2 (space) SNP Index y2 (space) \# of reads that indicate x2 and y2 belong to one community

For example, if you have 4 reads, the first read indicates SNP1 and SNP2 are in the same community. The second read indicates SNP1 and SNP2 are in different community. The third read indicate SNP1 and SNP3 are in the same community. The fourth read indicates SNP1 and SNP2 are in different community. Then the contact map should be:

3 3 3

1 1 -1     (The first read +1, the second read -1, the fourth read -1.)

1 3 1      (The third read +1.)

A demo contactmap file is provided in `contactmap_demo.txt`.


### Output File Format

We provide two formats of output. By default, the output is in blocks. Each block is printed separately. If you want a contiguous SNP sequence from start to the end ignoring the blocks, please use -nb or --noblock option when running.

### Algorithm procedure

When we read in a refhap reads file, we will get the contactmap from all the reads and the quality scores. Then the blocks are determined by finding all the connected components from the contact map. Then spectral-stitching algorithm runs on each block seperately. In the future release we will add parallel property for our program.


Note: In the new implementation of python version, we use the same definition of block as refhap and probhap etc, which is not identical to 10X definition, the one we use for evaluation package. For this version we are showing less blocks. The comparison of two different version is shown in `evaluation\new_implemention_comparison.png`.


## Usage of evaluation package

In folder `evaluation\`, we provide a demo for running on chromosome 20 of NA12878 WGS data.

First run `subsample&preprocess.py` for subsampling the coverage and generate adjacent matrix file required for main program.

Then run `pipeline.m` to run spectral-stitching algorithm on generated adjacent matrix. 

Finally, run `measure_performance_across_subsamples.R` to get the metrics switch error rate and unphased SNPs. Run `N50.py` to get N50.

If you have a a adjacent matrix (or contact map) in the same format as `contactmap_demo.csv`, You can run in MATLAB

```
phased_seq = Spectral_stitching(`contactmap_demo.csv`);
```

to get the phased sequence.
