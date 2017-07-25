# Spectral-Stitching for Haplotype Phasing

This repository implements the Spectral-Stitching algorithm for haplotype phasing proposed in:

* Yuxin Chen, Govinda Kamath, Changho Suh, and David Tse,  "[Community recovery in graphs with locality](http://proceedings.mlr.press/v48/chena16.html)," *International Conference on Machine Learning*, pp. 689-698, June 2016.
  * [Paper (ICML)](http://www.princeton.edu/~yc5/publications/Locality_ICML.pdf)
  * [Paper (Arxiv)](https://arxiv.org/abs/1602.03828)
  * [Slides](http://www.princeton.edu/~yc5/slides/Locality_ICML_slides.pdf)  

We provide a python main program, along with an evaluation program in MATLAB and R for running and evaluating on NA12878 WGS data in 10X paper.

Our algorithm is efficient and powerful especially when coverage is relatively deep (at least > 5X). When the coverage is shallow, our algorithm degrades into the simple spectral algorithm. 

The main program uses [NumPy](http://www.numpy.org/) and [SciPy](https://www.scipy.org/) packages, and accepts input and output in two formats. The evaluation code requires [CVX](http://cvxr.com/cvx/) in MATLAB as well as foreach & doParallel package in R. 


If you have any question, please feel free to contact 

[Banghua Zhu](mailto:13aeon.v01d@gmail.com), Tsinghua University   
[Yuxin Chen](mailto:yuxin.chen@princeton.edu), Princeton University  
[David Tse](mailto:dntse@stanford.edu), Stanford University


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






