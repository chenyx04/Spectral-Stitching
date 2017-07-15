Sys.time() # to have an idea of time passed

# set up parallel setting
items <- numSNP
detectedCores <- 4
batches <- detectedCores * 40
batchSets <- split(seq(1,numSNP), rep(1:batches, length.out=items))

cl<-makeCluster(4, outfile="parallelProgress")
registerDoParallel(cl)

# check for unphased individual SNPs
final_res <- foreach(p = seq(1,length(batchSets)),.combine = "comb",
                     .packages='Matrix',.inorder=FALSE) %dopar% {

    print(p)
    report_sdp <- NULL
    report_spec <- NULL

    for (i in seq(1,length(batchSets[[p]]))){

      v <- batchSets[[p]][i]   # actual index of the SNP

      # which reads cover this SNP
      ii = coveringReads[[v]]

      check = 0

      if (length(ii)!=0){      # no read covers this SNP
        check = 1
        count_sdp = 0; count_spec = 0
        # flip at that particular loci
        sdp_res_flip = sdp_res; spec_res_flip = spec_res;
        #if (sdp_res[v]==1){sdp_res_flip[v]=2} else {sdp_res_flip[v]=1}
        if (spec_res[v]==1){spec_res_flip[v]=2} else {spec_res_flip[v]=1}

        for (j in seq(1,length(ii))){  # loop over linakge reads covering it

          # remove SNPs that are too far away cause we never used them
          SNP_info <- SNP_covered[[ii[j]]]
          qq <- (abs(SNP_info-v)<=183)    # parameter to set
          SNP_info <- SNP_info[qq]
          res_info <- read_res[[ii[j]]]
          res_info <- res_info[qq]

          if (length(SNP_info) != 1){

            check = check+1   # assert that there's linkage read covering this SNP

            #sdp <- sdp_res[SNP_info]
            spec <- spec_res[SNP_info]
            #sdp_flip <- sdp_res_flip[SNP_info]
            spec_flip <- spec_res_flip[SNP_info]

            # does flipping help for each linkage read relationship
            # if (check_consistency(sdp,res_info) > check_consistency(sdp_flip,res_info)){
            #   count_sdp = count_sdp+1  # helps to flip
            # }else if (check_consistency(sdp,res_info) < check_consistency(sdp_flip,res_info)){
            #   count_sdp = count_sdp-1  # doesn't help to flip
            # }

            if (check_consistency(spec,res_info) > check_consistency(spec_flip,res_info)){
              count_spec = count_spec+1  # helps to flip
            }else if (check_consistency(spec,res_info) < check_consistency(spec_flip,res_info)){
              count_spec = count_spec-1  # doesn't help to flip
            }
          }

        }
      }

      # not sure if flipping and non-flipping at this SNP are equally likely
      # report contains unphased SNPs
      # they did it with abs(count_sdp) <= 1
      #if (check==0){report_sdp = c(report_sdp,v)}
      #if (check==0){report_spec = c(report_spec,v)}
      
      #if ( check==0){report_sdp = c(report_sdp,v)}
     if ( (check !=0 && abs(count_spec)<=1) || check==0){report_spec = c(report_spec,v)}

   }

    return(list(report_sdp,report_spec))
}
#report_sdp <- sort(final_res[[1]],decreasing=F)
report_spec <- sort(final_res[[2]],decreasing=F)

cat('Processing chr', chrNum, 'with coverage = ', coverage, '\n')

# percentage not phased
cat('Overall SNP num: ', numSNP, '\n')
#cat('PRELIMINARY: sdp % SNPs unphased ', length(report_sdp)/numSNP, '\n')
cat('PRELIMINARY: spec % SNPs unphased ', length(report_spec)/numSNP, '\n')
#print('unphased individual SNP done')


# set up parallel setting
items <- numSNP-1
detectedCores <- 4 #change to your number of cores
batches <- detectedCores * 40
batchSets <- split(seq(1,numSNP-1), rep(1:batches, length.out=items))

# determine breakpoints for phase blocks
final_res_1 <- foreach(p = seq(1,length(batchSets)),.combine = "comb",
                     .packages='Matrix',.inorder=FALSE) %dopar% {

     print(p)

     report_sdp_phase <- NULL
     report_spec_phase <- NULL

     for (i in seq(1,length(batchSets[[p]]))){

         v <- batchSets[[p]][i]   # actual index of the SNP
         check = 0; check1 = 0;
         count_sdp = 0; count_spec = 0;

         # remove the unphased SNPs from further consideration - SDP method
#          if (!(v %in% report_sdp)){

#            # only involves reads that connect left and right
#              ii <- which(sapply(SNP_covered, function(x){any(x>v) && any(x<=v)}) == TRUE)
# #ii = bridgingReads[[v]]

#            if (length(ii)!=0){      # no read covers this SNP

#              # flip everything to the right
#              left_sdp <- sdp_res[1:v]
#              right_sdp <- sdp_res[-v]
#              right_flip = right_sdp
#              right_flip[right_sdp==1] = 2; right_flip[right_sdp==2] = 1
#              sdp_res_flip = c(left_sdp, right_flip)


#              for (j in seq(1,length(ii))){

#                # remove SNPs that are too far away cause we never used them
#                # analogous to the "clustering step"
#                SNP_info <- SNP_covered[[ii[j]]]
#                qq <- (abs(SNP_info-v)<=183)   # parameter to set
#                SNP_info <- SNP_info[qq]
#                res_info <- read_res[[ii[j]]]
#                res_info <- res_info[qq]

#                # remove SNPs that we don't need to phase
#                pp <- which(SNP_info %in% report_sdp)
#                if (length(pp)!=0){
#                  SNP_info <- SNP_info[-pp]
#                  res_info <- res_info[-pp]
#                }

#                # does it still cover both sides? / not singleton read
#                if (length(SNP_info) != 1 && (any(SNP_info > v) && any(SNP_info<=v) == TRUE)){

#                  check = check+1

#                  sdp <- sdp_res[SNP_info]
#                  sdp_flip <- sdp_res_flip[SNP_info]

#                  # does flipping help for each linkage read relationship
#                  left <- check_consistency(sdp,res_info)
#                  right <- check_consistency(sdp_flip,res_info)
#                  if (left > right){
#                    count_sdp = count_sdp+(left-right)
#                  } else if (left < right) {
#                    count_sdp = count_sdp-(right-left)
#                  }

#                }

#              }
#            }

#            # not sure if flipping and non-flipping at this SNP are equally likely: report as unphased
#            if (check==0){
#              report_sdp_phase = c(report_sdp_phase,v)
#            }

# 	}

      # remove unphased SNPs from furthur consideration - spectral method
      if (!(v %in% report_spec)){

           # only involves reads that connect left and right
           ii <- which(sapply(SNP_covered, function(x){any(x>v) && any(x<=v)}) == TRUE)

           if (length(ii)!=0){      # no read covers this SNP

             # flip everything to the right
             left_spec <- spec_res[1:v]
             right_spec <- spec_res[-v]
             right_flip = right_spec
             right_flip[right_spec==1] = 2; right_flip[right_spec==2] = 1
             spec_res_flip = c(left_spec, right_flip)

             for (j in seq(1,length(ii))){

               # remove SNPs that are too far away cause we never used them
               # analogous to the "clustering step"
               SNP_info <- SNP_covered[[ii[j]]]
               qq <- (abs(SNP_info-v)<=183)   # parameter to set
               SNP_info <- SNP_info[qq]
               res_info <- read_res[[ii[j]]]
               res_info <- res_info[qq]

               # remove SNPs that we don't need to phase
               pp <- which(SNP_info %in% report_spec)
               if (length(pp)!=0){
                 SNP_info <- SNP_info[-pp]
                 res_info <- res_info[-pp]
               }

               # does it still cover both sides? / not singleton read
               if (length(SNP_info) != 1 && (any(SNP_info > v) && any(SNP_info<=v) == TRUE)){

                 check1 = check1+1

                 spec <- spec_res[SNP_info]
                 spec_flip <- spec_res_flip[SNP_info]

                 # does flipping help for each linkage read relationship
                 left <- check_consistency(spec,res_info)
                 right <- check_consistency(spec_flip,res_info)
                 if (left > right){
                   count_spec = count_spec + (left-right)
                 } else if (left < right) {
                   count_spec = count_spec - (right-left)
                 }

               }

             }
           }

           # not sure if flipping and non-flipping at this SNP are equally likely: report as unphased
           if ( (check1!=0 && abs(count_spec)<=1) || check1==0){
             report_spec_phase = c(report_spec_phase,v)
           }

	}

    }

    return(list(report_sdp_phase,report_spec_phase))
}

stopCluster(cl)


#report_sdp_phase = sort(final_res_1[[1]],decreasing=F)
report_spec_phase = sort(final_res_1[[2]],decreasing=F)
report_spec_phase = c(1,report_spec_phase)
report_spec_phase = c(report_spec_phase,length(spec_res))

# how many phase blocks are there
#print("Number of phase blocks: sdp, spec")
#numPhaseBlocks_sdp = length(report_sdp_phase-1)
numPhaseBlocks_spec = length(report_spec_phase-1)
#cat('numPhaseBlocks_sdp: ', numPhaseBlocks_sdp, '\n')
cat('Num of Phaseed Blocks for spec: ', numPhaseBlocks_spec, '\n')
# longest phase block
#print("Length of longest phase block: sdp, spec")
#longest_sdp_phase <- max(diff(pos[report_sdp_phase]));
longest_spec_phase <- max(diff(pos[report_spec_phase]));
#cat('longest_sdp_phase: ', longest_sdp_phase, '\n')
cat('Num of longest spec phase: ', longest_spec_phase, '\n')

#print('phase block breakpoint done')


# compute long and short switch error for phased SNPs

# check phase block length
#t1 <- which(diff(report_sdp_phase)==1)
#singleton_sdp <- report_sdp_phase[t1]
t2 <- which(diff(report_spec_phase)==1)
singleton_spec <- report_spec_phase[t2]

# singleton phase block is essentially an unphased SNP
#report_sdp_1 <- unique(c(report_sdp,singleton_sdp))
report_spec_1 <- unique(c(report_spec,singleton_spec))
# udpate the phasing percentage
#print('Updated % SNPs unphased: sdp, spec')
#cat('Updated sdp % SNPs unphased: ', length(report_sdp_1)/numSNP,'\n')
cat('Updated spec % SNPs unphased: ', length(report_spec_1)/numSNP,'\n')


# calculate error for each phase block
# err_sdp_short = 0; err_sdp_long = 0; len_r = 0; err_sdp_tradition = 0;

# for (j in seq(1,length(report_sdp_phase)-1)){
  
#   ll <- report_sdp_phase[j]
#   rr <- report_sdp_phase[j+1]-1
  
#   temp_res <- sdp_res[ll:rr]
#   temp_truth <- ground_truth[ll:rr]
  
#   if (length(temp_res) != 1){     # not singleton phase block
    
#     # remove unphased SNPs in the block
#     to_remove <- which(seq(ll,rr) %in% report_sdp)
#     if (length(to_remove)!=0){
#       temp_res <- temp_res[-to_remove]
#       temp_truth <- temp_truth[-to_remove]
#     }
    
#     if (length(temp_truth)!=1){   # after removing unphased SNP, not singleton
      
#       # flip to figure out the optimal global phase within phase block
#       our = temp_res
#       our[temp_res==1]=2
#       our[temp_res==2]=1
#       if (sum(temp_res!=temp_truth) > sum(our!=temp_truth)){
#         temp_res = our
#       }
#       mis_classify <- as.numeric(temp_res!=temp_truth)    
      
#       # find optimal decomposition into short & long switch error   
#       rle_res <- rle(mis_classify)
#       err_count <- rle_res$lengths[which(rle_res$values==1)]
#       err_sdp_tradition <- err_sdp_tradition + length(rle_res$lengths) - 1
#       len_r <- len_r+length(mis_classify)
#       err_sdp_short = err_sdp_short+sum(err_count[err_count<5])
#       err_sdp_long = err_sdp_long+length(which(err_count>=5)) 
#     }
    
#   }
# }


# sdp_short <- err_sdp_short/len_r
# sdp_long <- err_sdp_long/len_r
# sdp_tradition <- err_sdp_tradition/len_r
# cat("sdp len: ", len_r, "\n")
# cat("sdp short switch error rate: ", sdp_short, "\n")
# cat("sdp long switch error rate: ", sdp_long, "\n")
# cat("sdp tradition switch error rate: ", sdp_tradition, "\n")


# spectral method switch error
err_spec_short = 0; err_spec_long = 0; len_r = 0; err_spec_tradition = 0;

for (j in seq(1,length(report_spec_phase)-1)){
  
  ll <- report_spec_phase[j]
  rr <- report_spec_phase[j+1]-1
  
  temp_res <- spec_res[ll:rr]
  temp_truth <- ground_truth[ll:rr]
  
  if (length(temp_res) != 1){     # not singleton phase block
    
    # remove unphased SNPs in the block
    to_remove <- which(seq(ll,rr) %in% report_spec)
    if (length(to_remove)!=0){
      temp_res <- temp_res[-to_remove]
      temp_truth <- temp_truth[-to_remove]
    }
    
    if (length(temp_truth)!=1){   # after removing unphased SNP, not singleton
      
      # flip to figure out the optimal global phase within phase block
      our = temp_res
      our[temp_res==1]=2
      our[temp_res==2]=1
      if (sum(temp_res!=temp_truth) > sum(our!=temp_truth)){
        temp_res = our
      }
      mis_classify <- as.numeric(temp_res!=temp_truth)    
      
      # find optimal decomposition into short & long switch error   
      rle_res <- rle(mis_classify)
      err_count <- rle_res$lengths[which(rle_res$values==1)]
      err_spec_tradition <- err_spec_tradition + length(rle_res$lengths) - 1
      len_r <- len_r+length(mis_classify)
      err_spec_short = err_spec_short+sum(err_count[err_count<5])
      err_spec_long = err_spec_long+length(which(err_count>=5)) 
      #      print(err_count)

    }
    
  }
}

spec_short <- err_spec_short/len_r
spec_long <- err_spec_long/len_r
spec_tradition <- err_spec_tradition /len_r
cat("spec len: ", len_r, "\n")
cat("spec short switch error rate: ", spec_short, "\n")
cat("spec long switch error rate: ", spec_long, "\n")
cat("spec tradition switch error rate: ", spec_tradition, "\n")
#print('short/long switch error done'\n)

# write computed variables to file
cur_wd = dirname(sys.frame(1)$ofile) 
setwd(paste(cur_wd,"/data/metric/",sep=""))
dirName = paste("chr", chrNum, "_cov", coverage, sep = '')
dir.create(dirName)
setwd(dirName)
#write(report_sdp, file = "report_sdp", ncolumns = 1)
write(report_spec, file = "report_spec", ncolumns = 1)
write(singleton_spec, file = "singleton_spec", ncolumns = 1)
#write(singleton_sdp, file = "singleton_sdp", ncolumns = 1)
#write(report_sdp_phase, file = "report_sdp_phase", ncolumns = 1)
write(report_spec_phase, file = "report_spec_phase", ncolumns = 1)

Sys.time() # to have an idea of time passed