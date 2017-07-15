'''sub-samples and pre-processes fragment information
Note: everywhere in this code read should be replaced by fragment for correct semantics; the entire analysis revolves arounf fragments, not individual reads'''

from __future__ import print_function
from random import shuffle
for numchr in [20]:
    chrNum = str(numchr)
    groundTruthDir = "data/groundtruth/"
    preProcessedDir = "data/adjacent/"

    readDict = {}; resRead = {};
    # for each SNP record the fragments that touch it and the phase they assign to it (true / false) (the delimiter for the phase is a ',')
    # also create the reverse mapping: record for each fragment which SNPs it touches
    with open(groundTruthDir + 'chr' + chrNum + '_ground_truth.txt') as f:
        SNP_reads = f.readlines()
        numSNP = len(SNP_reads)
        for i, line in enumerate(SNP_reads):
            reads = line.split()[7]
           # print(i)
            if ',' in reads:
                # split fragments based on the phase they assign to the SNP
                true_reads, false_reads = reads.split(',')
                true_set, false_set = set((r.split('-')[0] for r in true_reads.split(';'))), set((r.split('-')[0] for r in false_reads.split(';')))
                unexpected_intersection = true_set & false_set
                true_set -= unexpected_intersection
                false_set -= unexpected_intersection
                assert len( true_set & false_set ) == 0
            else:
                true_reads = reads
                false_reads = None
                true_set = set((r.split('-')[0] for r in true_reads.split(';')))
            
            for trueRead in true_set:
                if trueRead in readDict:
                    readDict[trueRead].append(i+1) # 1-indexing
                    resRead[trueRead].append(True)
                else:
                    readDict[trueRead] = [i+1] # 1-indexing
                    resRead[trueRead] = [True]

            if false_reads == None:
                continue
            for falseRead in false_set:
                if falseRead in readDict:
                    readDict[falseRead].append(i+1) # 1-indexing
                    resRead[falseRead].append(False)
                else:
                    readDict[falseRead] = [i+1] # 1-indexing
                    resRead[falseRead] = [False]

    #print('num SNPs:', numSNP)
    #print('num Reads:', len(readDict))

    # remove fragments covering only one SNP (useless for phasing)
    toRemove = []
    for read in readDict:
        if len(readDict[read]) == 1:
            toRemove.append(read)
    for read in toRemove:
        del readDict[read]
        del resRead[read]
    numReads = len(readDict)
    #print('num Reads after deletion:', numReads)

    # sub-sample fragments according to the coverage ratio
    for coverage in [10, 13, 17, 23, 26, 37]:
        print('coverage:', coverage)
        subsampleReads = int(numReads * coverage / 37)

        # compute and store SNP_covered, res_read and coveringReads for the randomly selected sub-samples
        shuffleHelper = range(len(readDict))
        readDictVal = readDict.values()
        #print(shuffleHelper)
        shuffle(shuffleHelper)
        shuffleHelper = shuffleHelper[1:subsampleReads]
        #print(shuffleHelper)
        SNP_covered = []
        for i in shuffleHelper:
            SNP_covered.append(readDictVal[i])          # for every fragment i, SNP_covered[i] gives the ids of all the SNPs that this fragment touches

        resReadVal = resRead.values()
        read_res = []
        for i in shuffleHelper:
            read_res.append(resReadVal[i])              # for every fragment i, read_res[i] gives the phase (T/F) of all the SNPs that this fragment touches
            
        coveringReads = [[] for _ in range(numSNP)]     # for every SNP j, coveringReads gives the ids of the fragmnets that touch SNP j
        for i, read in enumerate(SNP_covered):
            for SNP in read:
                coveringReads[SNP-1].append(i+1)        # respect 1- and 0-indexing

        with open(preProcessedDir +  'chr' + chrNum + 'SNP_covered_cov' + str(coverage) + 'x', 'w') as f:
            for listSNP in SNP_covered:
                print( *listSNP, file = f)
                
        with open(preProcessedDir +  'chr' + chrNum + 'read_res_cov' + str(coverage) +'x', 'w') as f:
            for item in read_res:
                print( *item, file = f)
                
        with open(preProcessedDir +  'chr' + chrNum + 'coveringReads_cov' + str(coverage) +'x', 'w') as f:
            for item in coveringReads:
                print( *item, file = f)

        # for each SNP, store the bridging fragments: fragments that touch both a SNP at its left and a SNP at its right (not necessarily consecutive)
        bridgingReads = [[] for _ in range(numSNP)]
        for i, read in enumerate(SNP_covered):
            for SNPindex in range(1, len(read)-1):
                bridgingReads[read[SNPindex]-1].append(i+1) # respect 1- and 0-indexing
                
        with open(preProcessedDir +  'chr' + chrNum + 'bridgingReads_cov' + str(coverage) +'x', 'w') as f:
            for item in bridgingReads:
                print( *item, file = f)

        # compute and store adjoint matrices
        # for every fragment that links 2 SNPs i, j we add 1 to the corresponding entry of the adjoint matrix
        # (adj[i][j]) if thay share the same phase according to the fragment; otherwise we subtract 1
        adj = {}
        for readId, read in enumerate(SNP_covered):
            for i in range(len(read)):
                for j in range(i+1, len(read)):
                    assert read[i] < read[j]
                    pair = (read[i], read[j])
                    if pair in adj:
                        adj[pair] += 1 if read_res[readId][i] == read_res[readId][j] else -1
                    else:
                        adj[pair] = 1 if read_res[readId][i] == read_res[readId][j] else -1

        # store the matrix
        numLinks = len(adj)
        print('numLinks:', numLinks)
        with open(preProcessedDir +  'chr' + chrNum + 'adj_updated_cov' + str(coverage) + '.csv', 'w') as f:
            print("%%MatrixMarket matrix coordinate integer symmetric", file = f)
            print(numSNP, numSNP, numLinks, file = f)
            for i in range(numSNP):
                for j in range(i+1, numSNP):
                    # only store non-zero values of the adjoint matrix (it is very sparse)
                    if (i,j) in adj and adj[(i,j)] != 0:
                        print(j, i, adj[(i,j)], file = f)
