'''Computes the N50 metric across subsamples'''

from __future__ import print_function
chrNum = "20"

def blockLength(b, phase):
        '''Computes the length of phase block b'''
        i = phase.index(b)
        #print(b,i,phase[i+1])
        return pos[phase[i+1]-1] - pos[b]
def numSNPs(b, phase):
        '''Computes the number of SNPs in phase block b'''
        i = phase.index(b)
        return phase[i+1] - b

def N50(phase):
        '''Computes the N50 metric'''
        blocksSortedByLength = sorted(phase[:-1], key = lambda x: blockLength(x,phase), reverse = True)

        # keeps going through the phase blocks from longest to shortest and accumulates the number of SNPs so far
        s = 0
        for block in blocksSortedByLength:
                s += numSNPs(block, phase)
                if s >= halfSNP:
                        # as soon as we accumulate half the SNPs, we have arrived at the N50 block length
                        return blockLength(block,phase)

# read in the positions of each SNP
pos = []
posFile = "data/groundtruth/" + 'chr' + chrNum + '_ground_truth.txt'
with open(posFile) as f:
    for line in f.readlines():
        pos.append(int(line.split()[1]))

numSNP = len(pos)
halfSNP = numSNP / 2

# for all sub-sample sizes
for coverage in [ 10, 13,17, 23, 26, 37]:#
    dirName = "data/metric/chr" + chrNum + "_cov" + str(coverage) + "/"

    # read in the phase block information
    # sdpFile = "report_sdp_phase"
    # sdpPhase = []
    # with open(dirName + sdpFile) as f:
    #     for line in f.readlines():
    #         sdpPhase.append(int(line.strip()))
    
    specFile = "report_spec_phase"
    specPhase = []
    with open(dirName + specFile) as f:
        for line in f.readlines():
            specPhase.append(int(line.strip()))
        
        # compute the N50 metric
   # N50_sdp = N50(sdpPhase)
    N50_spec = N50(specPhase)

    print("coverage", coverage, "N50")
    print( '\tspec:', N50_spec)
