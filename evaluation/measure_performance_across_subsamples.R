rm(list = ls())
library('Matrix')
library('Matrix')
library('foreach')
library('doParallel')

# subfunction to check misclassification error
check_consistency <- function(our,res){
  
  res[res==TRUE] = 1; res[res==FALSE] = 2;
  # check how many SNP result doesn't agree with linkage read info
  temp = our 
  temp[our==1]=2
  temp[our==2]=1
  
  if (sum(our!=res) > sum(temp!=res)){
    return(sum(temp!=res))     # misclassification error
  } else{
    return(sum(our!=res))
  }
  
}

# define combine function
comb <- function(x, ...){  
  mapply(c,x,...,SIMPLIFY=FALSE)
}

chrNum = "20" # change this

for (coverage in c(10, 13, 17,23,26,37)) {#

  # load pre-processed linked read data
  cur_wd = dirname(sys.frame(1)$ofile) # Select your working directory
  setwd(paste(cur_wd, "/data/adjacent/",sep=""))
  SNP_covered = lapply(strsplit(readLines(
    con = paste("chr", chrNum, "SNP_covered_cov", coverage, 'x', sep = '')), " "), strtoi)
  coveringReads = lapply(strsplit(readLines(
    con = paste("chr", chrNum, "coveringReads_cov", coverage, 'x', sep = '')), " "), strtoi)
  read_res = lapply(strsplit(readLines(
    con = paste("chr", chrNum, "read_res_cov", coverage, 'x', sep = '')), " "), as.logical)
  bridgingReads = lapply(strsplit(readLines(
    con = paste("chr", chrNum, "bridgingReads_cov", coverage, 'x', sep = '')), " "), strtoi)
  
  setwd(paste(cur_wd,"/data/groundtruth/",sep=""))
  
  # load SNP number, their position and the ground truth phase
  ground_truth <- as.matrix(read.table(paste('chr', chrNum, '_ground_truth.txt', sep = ''),
                                       colClasses = c(rep("NULL",8), "integer")), mode="numeric")
  pos <- as.matrix(read.table(paste('chr', chrNum, '_ground_truth.txt', sep = ''),
                              colClasses = c("NULL", "integer", rep("NULL",7))), mode="numeric")
  numSNP = length(pos)
  
  # load our results
  setwd(paste(cur_wd,"/data/res/",sep=""))
  our_result <- read.table(paste('chr', chrNum, '_cov', coverage, sep = ''),sep=",")
  spec_res <- as.vector(our_result[1,], mode="numeric")

  
  # calculate metrics
  setwd(cur_wd)
  source("metrics_parallel.R",encoding='utf-8')
}
