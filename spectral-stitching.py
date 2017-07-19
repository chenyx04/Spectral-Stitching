#!/usr/bin/python
import argparse
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from math import exp, log
from lib.spectral import spectral
import time

##########              Argument Setting                 ####################
#############################################################################

parser = argparse.ArgumentParser(description="Spectral-stitching, haplotype assembly using community detection in graph with locality")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-r','--reads', help='path to reads file, this is mutually exclusive with contact map file')
group.add_argument('-c','--contactmap', help='path to contact map, this is mutually exclusive with reads file')
parser.add_argument('-o','--outphase', help='output file path', required=True)
parser.add_argument('-nb','--noblock', help='simply output all snps instead of blocks', action='store_true')
args = parser.parse_args()



##########                Parsing Input                  ####################
#############################################################################

start_time = time.clock()
if args.reads is not None:
    # Transfer refhap reads to contact map file
    with open(args.reads) as f:
        num_reads, num_snp = (int(x) for x in f.readline().split())
        # Build up contact map from given reads
        contact_map = sp.lil_matrix((num_snp,num_snp),dtype=np.float)
        contact_map_cnt = sp.lil_matrix((num_snp,num_snp),dtype=np.int)
        print ("Parsing Refhap file...")
        for i,line in enumerate(f):
            fields = line.strip().split()
            if len(fields) == 2: continue
            # extract data from the read:
            covered_positions = list()
            alleles = dict()
            qscores = dict()

            qscore_list = list(fields[-1][::-1])
            num_pieces = int(fields[0])

            # Iterate over contiguous SNP groups in a single read
            for k in xrange(num_pieces):
                start = int(fields[2+2*k])-1
                allele_list = fields[2+2*k+1]

                # Iterate over each SNP in a contiguous SNP grouop
                for j, allele in zip(xrange(start,start+len(allele_list)), allele_list):
                    alleles[j] = [int(allele)]
                    x = qscore_list.pop()
                    qscores[j] = float(ord(x)-33)
                    covered_positions.append(j)

            # Ignore reads that only cover a single SNP
            if len(covered_positions) == 0:
                continue

            # Get probability of two SNP being the same community
            # P(same community) = \sum_i (P(same community | read_i) * P(read_i))
            for k in alleles.keys():
                for j in alleles.keys():
                    if j > k:
                        if alleles[j]==alleles[k]:
                            contact_map[k,j] = contact_map[k,j] + (1 - 10**(-(qscores[k])/10.0)) * (1 - 10**(-(qscores[j])/10.0))
                        else:
                            contact_map[k,j] = contact_map[k,j] + (1 - (1 - 10**(-(qscores[k])/10.0)) * (1 - 10**(-(qscores[j])/10.0)))
                        contact_map_cnt[k,j] = contact_map_cnt[k,j] + 1
        # Transfer probability into 1 (same community) and -1 (different community)
        cm_row, cm_col = contact_map_cnt.nonzero()
        for idx in range(len(cm_row)):
            if contact_map[cm_row[idx],cm_col[idx]]/contact_map_cnt[cm_row[idx],cm_col[idx]]>0.5:
                contact_map[cm_row[idx],cm_col[idx]] = 1
            elif contact_map[cm_row[idx],cm_col[idx]]/contact_map_cnt[cm_row[idx],cm_col[idx]] == 0.5:
                contact_map[cm_row[idx],cm_col[idx]] = 0
            else:
                contact_map[cm_row[idx],cm_col[idx]] = -1
        
else:
    
    # Parse contact map file
    with open(args.contactmap) as f:
        print ("Parsing contactmap file...")
        num_snp, verbose, num_links = (int(x) for x in f.readline().split())        
        contact_map = sp.lil_matrix((num_snp,num_snp),dtype=np.float)
        if num_snp != verbose:
            print("Contact map size error! # of Rows != # of Columns!")
            exit(0)
        # Read in file
        for i,line in enumerate(f):
            fields = line.strip().split()
            if i == 0: continue     
            contact_map[int(fields[0])-1,int(fields[1])-1] = int(fields[2])
        # Transfer probability into 1 (same community) and -1 (different community)
        contact_map[contact_map>0]=1
        contact_map[contact_map<0]=-1   

##########                Preprocess Input               ####################
#############################################################################




# Make symmetric
contact_map_t = contact_map.transpose(copy=True)
contact_map = contact_map + contact_map_t

# Find all connected components
n_components, labels = connected_components(contact_map, directed=False, return_labels=True)
phased_snp = np.zeros((num_snp,))

out_file = open(args.outphase, "w")
out_file.write("")
out_file.close()

preprocess_time = time.clock()
print("Preprocessing finished! Elapsed time = " + str(preprocess_time - start_time) + 's' )


##########            Run Spectral-stitching             ####################
#############################################################################

W = 100
idx_dic = np.array(range(num_snp))
sub_idx_idc = np.array(range(W))
print ("Phasing " + str(n_components)+" blocks (including singleton SNP)...")
# Iterate over all connected components (blocks)
for i in range(n_components):

# Step1. Spectral
 
    # Extract position and elements of block i from contact map
    pos = idx_dic[labels==i]
    if len(pos) == 1: continue
    elements = contact_map[pos,:][:,pos]

    # Seperate elements into chunk with size W*W
    chunk_n = int(np.ceil((len(pos)-W)/W*2))+1
    # For chunk size <= 1, stitching is not needed
    if chunk_n <= 1:
        phased_snp[pos] = spectral(elements.toarray(), 2)
    else:
        chunk_snp={}
        # Iterate over all chunks
        for j in range(chunk_n):
            ww1 = j*W/2
            if j == chunk_n - 1:
                ww2 = len(pos)
            else:
                ww2 = j*W/2+W
            chunk = elements[ww1:ww2,:][:,ww1:ww2]
            # Seperate each chunk into blocks again, we run spectral only on connected components in chunk
            sub_n_components, sub_label = connected_components(chunk, directed=False, return_labels=True)            
            if sub_n_components>1:
                for sub_i in range(sub_n_components):
                    sub_pos = sub_idx_idc[sub_label==sub_i]
                    sub_chunk = chunk[sub_pos,:][:,sub_pos]
                    phased_snp[pos[ww1:ww2][sub_pos]] = spectral(sub_chunk.toarray(), 2)
            else:
                phased_snp[pos[ww1:ww2]] = spectral(chunk.toarray(), 2)
            chunk_snp[j] = phased_snp[pos[ww1:ww2]].flatten().copy()

        # Step2. Stitching

        for j in range(1, chunk_n):
            ww1 = j*W/2
            if j == chunk_n - 1:
                ww2 = len(pos)
            else:
                ww2 = j*W/2+W
            # Use the overlap region to determine flip or not
            if sum(chunk_snp[j-1][W/2:W] == chunk_snp[j][0:W/2]) < sum(chunk_snp[j-1][W/2:W] == -chunk_snp[j][0:W/2]):
                chunk_snp[j] = - chunk_snp[j]
            # For the overlap region, we determine the SNP by trusting the chunk with more reads included
            for k in range(W/2):
                if chunk_snp[j-1][W/2+k] != chunk_snp[j][k]:
                    if elements[ww1+k,:][:,(ww1-W/2):(ww1+W/2)].nnz >  elements[ww1+k,:][:,ww1:ww2].nnz:
                        chunk_snp[j][k] = chunk_snp[j-1][W/2+k]
            phased_snp[pos[ww1:ww2]] = chunk_snp[j] 

    # Step3. Local refinment

    for refine_iter in range(3):
        for tt, pos_tt in enumerate(pos):
            pp = elements[tt,:].nonzero()[1]
            pp = pp[abs(pp-tt)<(refine_iter+1)*W/2]
            info = elements[tt,:][:,pp].toarray().flatten()
            known_spec = phased_snp[pos[pp]].flatten()
            met = int(sum(info*known_spec))
            # Remain 0 if undetermined
            if met > 0:
                phased_snp[pos[tt]] = 1
            elif met < 0:
                phased_snp[pos[tt]] = -1

        ## Some immature thinking on long range switch error correction. Performance is not so good
        # for refine_iter in range(3):
        #     for tt, pos_tt in enumerate(pos):
        #         pp = elements[tt,:].nonzero()[1]
        #         pp = pp[abs(pp-tt)<(refine_iter+1)*W/2]
        #         left = pp[pp<tt]
        #         right = pp[pp>tt]
        #         info_left = elements[tt,:][:,left].toarray().flatten()
        #         known_spec_left = phased_snp[pos[left]].flatten()
        #         info_right = elements[tt,:][:,right].toarray().flatten()
        #         known_spec_right = phased_snp[pos[right]].flatten()
                
        #         met_left = int(sum(info_left*known_spec_left))
        #         met_right = int(sum(info_right*known_spec_right))
        #         # Remain 0 if undetermined
        #         if met_left > 5 and met_right < -5:
        #             phased_snp[pos[tt]] = 1
        #         elif met_left < -5 and met_right > 5:
        #             phased_snp[pos[tt]] = -1

spectral_time = time.clock()
print("Spectral-Stitching algorithm finished! Elapsed time = " + str(spectral_time - preprocess_time) +'s' )

##########                 Print Output                  ####################
#############################################################################

if args.noblock == True:
    # Do not print out blocks
    with open(args.outphase, "w") as testout:
        for i, snp in enumerate(phased_snp):
            if i == len(phased_snp)-1:
                if snp == 0:
                    snp = 1
                testout.write("%d" % (int(snp/2+1.5)))
            else:
                testout.write("%d," % (int(snp/2+1.5)))

else:
    # Print in the form of blocks
    out_file = open(args.outphase, "a")
    for i in range(n_components):
        pos = idx_dic[labels==i]
        if len(pos) == 1: continue
        out_file.write("BLOCK: offset: %d len: %d\n" % (pos[0]+1, len(pos)))
        for j in range(len(pos)):
            if phased_snp[pos[j]] == 0:
                out_file.write("%d \t -\n" % (pos[j]+1))
            else:
                out_file.write("%d \t %d\n" % (pos[j]+1, (1+phased_snp[pos[j]])/2))
        out_file.write("********\n")
    out_file.close()
