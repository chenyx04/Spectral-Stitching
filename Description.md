---
layout: default
title: Home
---
[Home](https://chenyx04.github.io/Spectral-Stitching/)  
[Description](https://chenyx04.github.io/Spectral-Stitching/Description)  
[Quick Start](https://chenyx04.github.io/Spectral-Stitching/users_guide)

# Haplotype Phasing


Humans have 23 pairs of homologous chromosomes, which are identical except on certain positions called single nucleotide polymorphisms (SNPs). A haplotype of an individual is the pair of sequences of SNPs on the two homologous chromosomes. Knowing the haplotypes of individuals can lead to a better understanding of the interplay of genetic variation and disease as well as better inference of human demographic history. The haplotype phasing problem is that of inferring haplotypes of individuals from high-throughput sequencing data. 

The Spectral-Stitching algorithm consists of three stages:

**Stage 1: node splitting and spectral estimation.**  Split all nodes into several overlapping subsets Vl (l ≥
1) of size W, such that any two adjacent subsets share W/2 common vertices. We choose the size W
of each Vl to be r for rings / lines, and on the order of davg for other graphs. We then run spectral
methods separately on each subgraph Gl
induced by Vl

{% raw %}
  $$a^2 + b^2 = c^2$$ --> note that all equations between these tags will not need escaping! 
 {% endraw %}









