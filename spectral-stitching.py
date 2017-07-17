#!/usr/bin/python
import argparse
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from math import exp, log
from lib.spectral import spectral
import time

#############################################################################

parser = argparse.ArgumentParser(description="Spectral-stitching, haplotype assembly using community detection in graph with locality")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-r','--reads', help='path to reads file, this is mutually exclusive with contact map file')
group.add_argument('-c','--contactmap', help='path to contact map, this is mutually exclusive with reads file')
parser.add_argument('-o','--outphase', help='output file path', required=True)
parser.add_argument('-nb','--noblock', help='simply output all snps instead of blocks', action='store_true')
args = parser.parse_args()

#############################################################################
start_time = time.clock()
# figure out the number of positions to phase:
if args.reads is not None:
       
    with open(args.reads) as f:
        num_reads, num_snp = (int(x) for x in f.readline().split())
        data_start, data_end = 0, num_snp - 1
        contact_map = sp.lil_matrix((num_snp,num_snp),dtype=np.float)
        print ("Parsing Refhap file...")
        for i,line in enumerate(f):
            fields = line.strip().split()
            if len(fields) == 2: continue
            # extract data for the read:
            covered_positions = list()
            alleles = dict()
            qscores = dict()

            qscore_list = list(fields[-1][::-1])
            num_pieces = int(fields[0])

            for k in xrange(num_pieces):
                start = int(fields[2+2*k])-1
                allele_list = fields[2+2*k+1]

                for j, allele in zip(xrange(start,start+len(allele_list)), allele_list):
                    alleles[j] = [int(allele)]
                    x = qscore_list.pop()
                    qscores[j] = float(ord(x)-33)
                    covered_positions.append(j)

            if len(covered_positions) == 0:
                continue

            for k in alleles.keys():
                for j in alleles.keys():
                    if j > k:
                        if alleles[j]==alleles[k]:
                            contact_map[k,j] = contact_map[k,j] + 10**(-(qscores[k]+qscores[j])/10.0)
                        else:
                            contact_map[k,j] = contact_map[k,j] - 10**(-(qscores[k]+qscores[j])/10.0)

else:
    
    with open(args.contactmap) as f:
        print ("Parsing contactmap file...")
        num_snp, verbose, num_links = (int(x) for x in f.readline().split())        
        contact_map = sp.lil_matrix((num_snp,num_snp),dtype=np.float)
        if num_snp != verbose:
            print("Contact map size error! # of Rows != # of Columns!")
            exit(0)
        for i,line in enumerate(f):
            fields = line.strip().split()
            if i == 0: continue     
            contact_map[int(fields[0])-1,int(fields[1])-1] = int(fields[2])

contact_map[contact_map>0]=1
contact_map[contact_map<0]=-1   
contact_map_t = contact_map.transpose(copy=True)
contact_map = contact_map + contact_map_t
n_components, labels = connected_components(contact_map, directed=False, return_labels=True)
phased_snp = np.zeros((num_snp,))

out_file = open(args.outphase, "w")
out_file.write("")
out_file.close()

preprocess_time = time.clock()
print("Preprocessing finished! Elapsed time = " + str(preprocess_time - start_time) + 's' )


# phase each locally phaseable block in turn

idx_dic = np.array(range(num_snp))
W = 100
sub_idx_idc = np.array(range(W))
print ("Phasing " + str(n_components)+" blocks...")
for i in range(n_components):
    # Step1. Spectral
    pos = idx_dic[labels==i]
    if len(pos) == 1: continue
    elements = contact_map[pos,:][:,pos]
    chunk_n = int(np.ceil((len(pos)-W)/W*2))
    if chunk_n <= 1:
        phased_snp[pos] = spectral(elements.toarray(), 2)
    else:
        chunk_snp={}
        for j in range(chunk_n):
            ww1 = j*W/2
            if j == chunk_n - 1:
                ww2 = len(pos)
            else:
                ww2 = j*W/2+W
            chunk = elements[ww1:ww2,:][:,ww1:ww2]
            sub_n_components, sub_label = connected_components(chunk, directed=False, return_labels=True)
            if sub_n_components>1:
                for sub_i in range(sub_n_components):
                    sub_pos = sub_idx_idc[sub_label==sub_i]
                    sub_chunk = chunk[sub_pos,:][:,sub_pos]
                    #print(sub_chunk.toarray())
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
            if sum(chunk_snp[j-1][W/2:W] == chunk_snp[j][0:W/2]) < sum(chunk_snp[j-1][W/2:W] == -chunk_snp[j][0:W/2]):
                chunk_snp[j] = - chunk_snp[j]
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
                met = int(sum(info*known_spec) > 0)
                if met > 0:
                    phased_snp[pos[tt]] = 1
                elif met < 0:
                    phased_snp[pos[tt]] = -1

spectral_time = time.clock()
print("Spectral-Stitching algorithm finished! Elapsed time = " + str(spectral_time - preprocess_time) +'s' )
if args.noblock == True:
    with open(args.outphase, "w") as testout:
        for i, snp in enumerate(phased_snp):
            if i == len(phased_snp)-1:
                if snp == 0:
                    snp = 1
                testout.write("%d" % (int(snp/2+1.5)))
            else:
                testout.write("%d," % (int(snp/2+1.5)))

else:
    out_file = open(args.outphase, "a")
    for i in range(n_components):
        # Step1. Spectral
        pos = idx_dic[labels==i]
        if len(pos) == 1: continue
        out_file.write("BLOCK: offset: %d len: %d\n" % (pos[0], len(pos)))
        for j in range(len(pos)):
            if phased_snp[pos[j]] == 0:
                out_file.write("%d \t -\n" % (pos[j]))
            else:
                out_file.write("%d \t %d\n" % (pos[j], phased_snp[j]))
        out_file.write("********\n")
    out_file.close()