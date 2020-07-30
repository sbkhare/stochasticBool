# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 15:24:08 2020

@author: Sikander
"""

import probabilisticBoolean as pb
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import pickle as pkl

N = 3 #sys.argv[1]
K = 2

def stg_thresholding(lutnum, incr):
    err_lst = np.arange(0.25, 0.5, incr)
    
    pbn = pb.PBNetwork(N, K) #Number of possible combinations of deterministic LUTs and adjacency matrixes: (N*math.factorial(K))*(math.factorial(N-1)/(math.factorial(K)*math.factorial(N-1-K)))**N
    pbn.createNetwork()
    pbn.errLUT(lutnum, 0)
    det_tup = pbn.deterministicLUT()
    stg_det = pbn.detSTG(show=True)
    plt.savefig("thresholding/N{0}_K{1}_detstg_lutnum{2}.png".format(N, K, lutnum))
    
    for err in err_lst:
        pbn.errLUT(lutnum, err)
        pmax, pdet, kr, ke, ks = pbn.deterministicLUT()
        stg_prob, p_thresh = pbn.stateTransitionGraph(pmax=pmax, thresh=(0.5)**N, show=True, save=False)
        plt.savefig("thresholding/N{0}_K{1}_edgprob_lutnum{2}_err{3}.png".format(N, K, lutnum, err))
        stg_thresh = pb.absoluteThreshold(stg_prob, thresh=(0.5)**N, show=True)
        plt.savefig("thresholding/N{0}_K{1}_stochstg_lutnum{2}_err{3}.png".format(N, K, lutnum, err))

if __name__=='__main__':
    stg_thresholding(10, 0.01)