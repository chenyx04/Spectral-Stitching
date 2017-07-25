---
layout: default
title: Description
---

[Home](https://chenyx04.github.io/Spectral-Stitching/)  
[Description](https://chenyx04.github.io/Spectral-Stitching/Description)  
[Quick Start](https://chenyx04.github.io/Spectral-Stitching/users_guide)

# Haplotype Phasing


Humans have 23 pairs of homologous chromosomes, one maternal and one paternal. Each pair are
identical sequences of nucleotides A,G,C,T’s except on certain documented positions called single nucleotide
polymorphisms (SNPs), or genetic variants. At each of these positions, one of the chromosomes takes on
one of A,G,C or T which is the same as the majority of the population (called the major allele), while the
other chromosome takes on a variant (also called minor allele). A haplotype of an individual is the pair of sequences of SNPs on the two homologous chromosomes. Knowing the haplotypes of individuals can lead to a better understanding of the interplay of genetic variation and disease as well as better inference of human demographic history. The advent of next generation sequencing technologies allows haplotype phasing by providing linking reads between multiple SNP locations. The haplotype phasing problem aims at inferring haplotypes of individuals from high-throughput sequencing data. 

One can formulate the problem of haplotype phasing as recovery of two communities of SNP locations,
those with the variant on the maternal chromosome and those with the variant on the paternal chromosome. Each pair of linking reads gives a noisy measurement of whether two (or more) SNPs have the variant on
the same chromosome or different chromosomes.


# Spectral-Stitching Algorithm

The Spectral-Stitching algorithm consists of three stages:

**Stage 1: node splitting and spectral estimation.**  Split all nodes into several *overlapping* subsets, and run spectral
methods separately on each subgraph induced by each vertex subset, in the hope of achieving approximate estimates for each subgraph. 

<br>
 
![Image of Spectral Stitching1](Stage1.png)

<br>

**Stage 2: stiching the estimates.**  The aim of this stage is to stitch together the outputs of Stage 1
computed in isolation for the collection of overlapping subgraphs, so as to ensure that they have matching global phases. 

<br>
 
![Image of Spectral Stitching2](Stage2.png)

<br>

**Stage 3: successive local refinement.**  Clean up all estimates using both backward and
forward samples in order to maximize recovery accuracy. This is achieved by running local majority
voting from the neighbors of each vertex until convergence. 

<br>

![Image of Spectral Stitching3](Stage3.png)

<br>

Details can be found in the following paper:

* Yuxin Chen, Govinda Kamath, Changho Suh, and David Tse,  "[Community recovery in graphs with locality](http://proceedings.mlr.press/v48/chena16.html)," *International Conference on Machine Learning*, pp. 689-698, June 2016.





