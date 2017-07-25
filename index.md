---
title: 'Spectral-Stitching algorithm'
---

# Spectral-Stitching for Haplotype Phasing

This repository implements the Spectral-Stitching algorithm for haplotype phasing proposed in:

* Yuxin Chen, Govinda Kamath, Changho Suh, and David Tse,  "[Community recovery in graphs with locality](http://proceedings.mlr.press/v48/chena16.html)," *International Conference on Machine Learning*, pp. 689-698, June 2016.
  * [Paper (ICML)](http://www.princeton.edu/~yc5/publications/Locality_ICML.pdf)
  * [Paper (Arxiv)](https://arxiv.org/abs/1602.03828)
  * [Slides](http://www.princeton.edu/~yc5/slides/Locality_ICML_slides.pdf)  
  * [Github Repo](https://github.com/chenyx04/Spectral-Stitching)

The main program is implemented in Python. An evaluation program in MATLAB and R is also provided for running and evaluating on NA12878 WGS data in 10X paper.

The Spectral-Stitching algorithm is efficient especially when the coverage is relatively deep (at least > 5X). When the coverage is shallow, our algorithm degrades into the simple spectral algorithm. 

The main program uses [NumPy](http://www.numpy.org/) and [SciPy](https://www.scipy.org/) packages, and accepts input and output in two formats. The evaluation code requires [CVX](http://cvxr.com/cvx/) in MATLAB as well as foreach & doParallel package in R. 


If you have any question, please feel free to contact 

[Banghua Zhu](mailto:13aeon.v01d@gmail.com), Tsinghua University   
[Yuxin Chen](mailto:yuxin.chen@princeton.edu), Princeton University  
[David Tse](mailto:dntse@stanford.edu), Stanford University








