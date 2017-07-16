#!/usr/bin/python

import argparse
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from math import exp, log
#############################################################################

parser = argparse.ArgumentParser(description="Spectral-stitching, haplotype assembly using community detection in graph with locality")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-r','--reads', help='path to reads file, this is mutually exclusive with contact map file')
group.add_argument('-c','--contactmap', help='path to contact map, this is mutually exclusive with reads file')
parser.add_argument('-p','--phase', help='output file path', required=True)

args = parser.parse_args()

#############################################################################

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
phased_snp = np.zeros((num_snp,1))

out_file = open(args.phase, "w")
out_file.write("")
out_file.close()

print ("Phasing blocks...")
print (n_components)
idx_dic = np.array(range(num_snp))
for i in range(n_components):
    #print(labels==i)
    pos = idx_dic[labels==i]
    elements = contact_map[pos,:][:,pos]
# phase each locally phaseable block in turn

N_blocks = len(hmm.blocks)
for i, block in enumerate(hmm.blocks):
    if i % 100 == 0:
      print ("Block %d/%d" % (i, N_blocks))
    start,end = block[0], block[-1]

    out_file = open(args.phase, "a")
    out_file.write("BLOCK: offset: %d len: %d positions: %d\n" % (start+1, end-start+1, len(block)))
    # print "BLOCK: offset: %d len: %d positions: %d" % (start+1, end-start+1, len(block))

    assignments = open(args.assignments, 'a')

    haplotype = hmm.run_viterbi(block)
    hmm.first_y = haplotype[0][0]

    right_segments = set()
    for (y,s) in haplotype: right_segments |= s

    left_clouds = set([M[i].name for j in block for i in hmm.segments_at_position[j]
                       if i not in right_segments and M[i][j] != -1])
    right_clouds = set([M[i].name for j in block for i in hmm.segments_at_position[j]
                        if i in right_segments and M[i][j] != -1])
    left_cloud_names = ','.join([cloud_name for cloud_name in left_clouds])
    right_cloud_names = ','.join([cloud_name for cloud_name in right_clouds])

    for c in left_clouds:
        assignments.write('%s\t0\n' % (c))

    for c in right_clouds:
        assignments.write('%s\t1\n' % (c))

    hmm.run_forwards(block)
    hmm.run_backwards(block)

    for k,(y,s) in enumerate(haplotype):
        j = block[k]
        out_file.write("%d\t" % (j))
        if j in hmm.uncovered_positions:
            out_file.write("-\t-\t")
        else:
            out_file.write("%d\t%d\t" % y)
        if k == 0:
            score = 0.5
        else:
            score = hmm.get_posterior_transition_prob_y_only(
                k, block, y, y_prev)
        emission_probability = exp(hmm.log_emission_prob(j, y, s))
        posterior_code = exp(hmm.get_posterior_code_y_only(k,block,y))
        out_file.write("%f\t%f\t%f\n" % (score,posterior_code,emission_probability))

        y_prev = y

    out_file.write("********\n")

    out_file.close()
    assignments.close()
