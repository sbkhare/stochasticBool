# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 13:39:21 2020

@author: Sikander
"""

import probabilisticBoolean as pb
import sys
import matplotlib.pyplot as plt
import numpy as np
import math

N = 3 #sys.argv[1]
K = 2 #sys.argv[2]

FULL_ENSEMBLE = 1
DET_BIASED = 1
PROB_BIASED = 1 
UNIFORM = 0
########################
ERR = 0

p_max_lst = []
p_det_lst = []
kr_lst= []
ke_lst= []
ks_lst = []
classes = []
kr_dets = []
ke_dets = []
ks_dets = []
dev1_lst = []
dev2_lst = []

if FULL_ENSEMBLE ==1:
    if DET_BIASED == 1:
        for i in range(5000):
            pbn = pb.PBNetwork(N, K) #Number of possible combinations of deterministic LUTs and adjacency matrixes: (2**(2**K))*(math.factorial(N-1)/(math.factorial(K)*math.factorial(N-1-K)))**N
            pbn.createNetwork()
            pbn.biasedLUT(det=True)
            p_max, p_det, kr_det, ke_det, ks_det = pbn.deterministicLUT()
            kr_dets.append(kr_det)
            ke_dets.append(ke_det)
            ks_dets.append(ks_det)
            p_max_lst.append(p_max)
            p_det_lst.append(p_det)
            kr, ke, ks, clas, lut_distr = pbn.canalization(p_max)
            classes.append(clas)
            kr_lst.append(kr)
            ke_lst.append(ke)
            ks_lst.append(ks)
            stg_prob, p_thresh = pbn.stateTransitionGraph(p_max)
            stg_det = pbn.detSTG()
            dev1, dev2 = pb.stg_deviation(N, stg_det, stg_prob)
            dev1_lst.append(dev1)
            dev2_lst.append(dev2)
        #    stg_thresh = pb.cumulativeThreshold(stg_prob, thresh=0.5)
    
    if PROB_BIASED == 1: 
        for i in range(5000):
            pbn = pb.PBNetwork(N, K) #Number of possible combinations of deterministic LUTs and adjacency matrixes: (2**(2**K))*(math.factorial(N-1)/(math.factorial(K)*math.factorial(N-1-K)))**N
            pbn.createNetwork()
            pbn.biasedLUT(det=False)
            p_max, p_det, kr_det, ke_det, ks_det = pbn.deterministicLUT()
            kr_dets.append(kr_det)
            ke_dets.append(ke_det)
            ks_dets.append(ks_det)
            p_max_lst.append(p_max)
            p_det_lst.append(p_det)
            kr, ke, ks, clas, lut_distr = pbn.canalization(p_max)
            classes.append(clas)
            kr_lst.append(kr)
            ke_lst.append(ke)
            ks_lst.append(ks)
            stg_prob, p_thresh = pbn.stateTransitionGraph(p_max)
            stg_det = pbn.detSTG()
            dev1, dev2 = pb.stg_deviation(N, stg_det, stg_prob)
            dev1_lst.append(dev1)
            dev2_lst.append(dev2)
        #    stg_thresh = pb.cumulativeThreshold(stg_prob, thresh=0.5)
        
    if UNIFORM == 1:
        incr = 0.1
        bins = np.arange(0.5, 1, incr)
        for i, b in enumerate(bins):
            if i == 0:
                continue
            elif i == len(bins)-1:
                for j in range(1000):
                    pbn = pb.PBNetwork(N, K) #Number of possible combinations of deterministic LUTs and adjacency matrixes: (2**(2**K))*(math.factorial(N-1)/(math.factorial(K)*math.factorial(N-1-K)))**N
                    pbn.createNetwork()
                    pbn.biasedRangeLUT(mn=bins[i-1], mx=b)
                    p_max, p_det, kr_det, ke_det, ks_det = pbn.deterministicLUT()
                    kr_dets.append(kr_det)
                    ke_dets.append(ke_det)
                    ks_dets.append(ks_det)
                    p_max_lst.append(p_max)
                    p_det_lst.append(p_det)
                    kr, ke, ks, clas, lut_distr = pbn.canalization(p_max)
                    classes.append(clas)
                    kr_lst.append(kr)
                    ke_lst.append(ke)
                    ks_lst.append(ks)
                    stg_prob, p_thresh = pbn.stateTransitionGraph(p_max)
                    stg_det = pbn.detSTG()
                    dev1, dev2 = pb.stg_deviation(N, stg_det, stg_prob)
                    dev1_lst.append(dev1)
                    dev2_lst.append(dev2)
            else:
                for j in range(1000):
                    pbn = pb.PBNetwork(N, K) #Number of possible combinations of deterministic LUTs and adjacency matrixes: (N*math.factorial(K))*(2**(2**K))*(math.factorial(N-1)/(math.factorial(K)*math.factorial(N-1-K)))**N
                    pbn.createNetwork()
                    pbn.biasedRangeLUT(mn=bins[i-1], mx=b)
                    p_max, p_det, kr_det, ke_det, ks_det = pbn.deterministicLUT()
                    kr_dets.append(kr_det)
                    ke_dets.append(ke_det)
                    ks_dets.append(ks_det)
                    p_max_lst.append(p_max)
                    p_det_lst.append(p_det)
                    kr, ke, ks, clas, lut_distr = pbn.canalization(p_max)
                    classes.append(clas)
                    kr_lst.append(kr)
                    ke_lst.append(ke)
                    ks_lst.append(ks)
                    stg_prob, p_thresh = pbn.stateTransitionGraph(p_max)
                    stg_det = pbn.detSTG()
                    dev1, dev2 = pb.stg_deviation(N, stg_det, stg_prob)
                    dev1_lst.append(dev1)
                    dev2_lst.append(dev2)
                for j in range(1000):
                    pbn = pb.PBNetwork(N, K) #Number of possible combinations of deterministic LUTs and adjacency matrixes: (N*math.factorial(K))*(2**(2**K))*(math.factorial(N-1)/(math.factorial(K)*math.factorial(N-1-K)))**N
                    pbn.createNetwork()
                    pbn.biasedRangeLUT(mn=bins[i-1]+incr/2, mx=b+incr/2)
                    p_max, p_det, kr_det, ke_det, ks_det = pbn.deterministicLUT()
                    kr_dets.append(kr_det)
                    ke_dets.append(ke_det)
                    ks_dets.append(ks_det)
                    p_max_lst.append(p_max)
                    p_det_lst.append(p_det)
                    kr, ke, ks, clas, lut_distr = pbn.canalization(p_max)
                    classes.append(clas)
                    kr_lst.append(kr)
                    ke_lst.append(ke)
                    ks_lst.append(ks)
                    stg_prob, p_thresh = pbn.stateTransitionGraph(p_max)
                    stg_det = pbn.detSTG()
                    dev1, dev2 = pb.stg_deviation(N, stg_det, stg_prob)
                    dev1_lst.append(dev1)
                    dev2_lst.append(dev2)
        #    stg_thresh = pb.cumulativeThreshold(stg_prob, thresh=0.5)
      
    plt.figure()
    plt.title("Stochasticity and State Transition Graphs, N={0} K={1}".format(N, K))
    plt.xlabel("Average p_max")
    plt.ylabel("Deviation between stochastic and deterministic STG edges")
    plt.scatter(p_max_lst, dev1_lst, color='blue')
    plt.savefig("full_ensemble/N{0}_K{1}_stgdev1_pmax.png".format(N, K))
    plt.figure()
    plt.xlabel("p_det")
    plt.ylabel("Deviation between stochastic and deterministic STG edges")
    plt.scatter(p_det_lst, dev1_lst, color='orange')
    plt.title("Stochasticity and State Transition Graphs, N={0} K={1}".format(N, K))
    plt.savefig("full_ensemble/N{0}_K{1}_stgdev1_pdet.png".format(N, K))
    
    plt.figure()
    plt.title("Stochasticity and State Transition Graphs, N={0} K={1}".format(N, K))
    plt.xlabel("Average p_max")
    plt.ylabel("Deviation between stochastic and deterministic STG non-edges")
    plt.scatter(p_max_lst, dev2_lst, color='blue')
    plt.savefig("full_ensemble/N{0}_K{1}_stgdev2_pmax.png".format(N, K))
    plt.figure()
    plt.xlabel("p_det")
    plt.ylabel("Deviation between stochastic and deterministic STG non-edges")
    plt.scatter(p_det_lst, dev2_lst, color='orange')
    plt.title("Stochasticity and State Transition Graphs, N={0} K={1}".format(N, K))
    plt.savefig("full_ensemble/N{0}_K{1}_stgdev2_pdet.png".format(N, K))
    
    
    plt.figure()
    plt.title("Stochasticity and Input Redundancy")
    plt.xlabel("Average p_max")
    plt.ylabel("Deviation between stochastic and deterministic k_r")
    plt.scatter(p_max_lst, np.array(kr_lst)-np.array(kr_dets), color='blue')
    plt.savefig("full_ensemble/K{0}_krdev_pmax.png".format(K))
    plt.figure()
    plt.xlabel("p_det")
    plt.ylabel("Deviation between stochastic and deterministic k_r")
    plt.scatter(p_det_lst, np.array(kr_lst)-np.array(kr_dets), color='orange')
    plt.title("Stochasticity and Input Redundancy")
    plt.savefig("full_ensemble/K{0}_krdev_pdet.png".format(K))
    
    plt.figure()
    plt.title("Stochasticity and Effective Connectivity")
    plt.xlabel("Average p_max")
    plt.ylabel("Deviation between stochastic and deterministic k_e")
    plt.scatter(p_max_lst, np.array(ke_lst)-np.array(ke_dets), color='blue')
    plt.savefig("full_ensemble/K{0}_kedev_pmax.png".format(K))
    plt.figure()
    plt.xlabel("p_det")
    plt.ylabel("Deviation between stochastic and deterministic k_e")
    plt.scatter(p_det_lst, np.array(ke_lst)-np.array(ke_dets), color='orange')
    plt.title("Stochasticity and Effective Connectivity")
    plt.savefig("full_ensemble/K{0}_kedev_pdet.png".format(K))
    
    plt.figure()
    plt.title("Stochasticity and Input Symmetry")
    plt.xlabel("Average p_max")
    plt.ylabel("Deviation between stochastic and deterministic k_s")
    plt.scatter(p_max_lst, np.array(ks_lst)-np.array(ks_dets), color='blue')
    plt.savefig("full_ensemble/K{0}_ksdev_pmax.png".format(K))
    plt.figure()
    plt.xlabel("p_det")
    plt.ylabel("Deviation between stochastic and deterministic k_s")
    plt.scatter(p_det_lst, np.array(ks_lst)-np.array(ks_dets), color='orange')
    plt.title("Stochasticity and Input Symmetry")
    plt.savefig("full_ensemble/K{0}_ksdev_pdet.png".format(K))
    
    plt.figure()
    plt.subplot(121)
    plt.hist(p_max_lst, bins=50, color='blue')
    plt.xlabel("Average p_max")
    plt.ylabel("Count")
    plt.subplot(122)
    plt.hist(p_det_lst, bins=100, color='orange')
    plt.xlabel("p_det")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("full_ensemble/K{0}_distr_pmaxdet.png".format(K))
    
    plt.figure()
    #Color and label quadrants for classes
    plt.scatter(kr_lst, ks_lst, color='g')
    plt.axvline(x=K/2, linestyle='-', color='r')
    plt.axhline(y=K/2, linestyle='-', color='r')
    plt.xlabel("k_r")
    plt.ylabel("k_s")
    plt.title("Input Redundancy vs. Input Symmetry")
    plt.tight_layout()
    plt.savefig("full_ensemble/K{0}_distr_kr_ks.png".format(K))
    
if ERR == 1:
    lutnum = 6
    lut = 'XOR'
    incr = 0.01
    err_lst = np.arange(0, 0.5, incr)
    rng = 10*(N*math.factorial(K))*(math.factorial(N-1)/(math.factorial(K)*math.factorial(N-1-K)))**N
    err_lst2 = []
    for error in err_lst:
        for i in range(int(rng)):
            pbn = pb.PBNetwork(N, K) #Number of possible combinations of deterministic LUTs and adjacency matrixes: (N*math.factorial(K))*(math.factorial(N-1)/(math.factorial(K)*math.factorial(N-1-K)))**N
            pbn.createNetwork()
            pbn.errLUT(lutnum, error)
            p_max, p_det, kr_det, ke_det, ks_det = pbn.deterministicLUT()
            kr_dets.append(kr_det)
            ke_dets.append(ke_det)
            ks_dets.append(ks_det)
            p_max_lst.append(p_max)
            p_det_lst.append(p_det)
            kr, ke, ks, clas, lut_distr = pbn.canalization(p_max)
            classes.append(clas)
            kr_lst.append(kr)
            ke_lst.append(ke)
            ks_lst.append(ks)
            stg_prob, p_thresh = pbn.stateTransitionGraph(p_max)
            stg_det = pbn.detSTG()
            dev1, dev2 = pb.stg_deviation(N, stg_det, stg_prob)
            dev1_lst.append(dev1)
            dev2_lst.append(dev2)
            err_lst2.append(error)
            
    plt.figure()
    plt.title("Stochasticity and State Transition Graphs, N={0} K={1} LUT={2}".format(N, K, lut))
    plt.xlabel("Error rate")
    plt.ylabel("Deviation between stochastic and deterministic STG edges")
    plt.scatter(err_lst2, dev1_lst, color='blue')
    plt.savefig("error_rate/N{0}_K{1}_stgdev1_lutnum{2}.png".format(N, K, lutnum))
    
    plt.figure()
    plt.title("Stochasticity and State Transition Graphs, N={0} K={1} LUT={2}".format(N, K, lut))
    plt.xlabel("Error rate")
    plt.ylabel("Deviation between stochastic and deterministic STG non-edges")
    plt.scatter(err_lst2, dev2_lst, color='orange')
    plt.savefig("error_rate/N{0}_K{1}_stgdev2_lutnum{2}.png".format(N, K, lutnum))
    
    
    plt.figure()
    plt.title("Stochasticity and Input Redundancy, LUT={0}".format(lut))
    plt.xlabel("Error rate")
    plt.ylabel("Deviation between stochastic and deterministic k_r")
    plt.scatter(err_lst2, np.array(kr_lst)-np.array(kr_dets), color='red')
    plt.savefig("error_rate/K{0}_krdev_lutnum{1}.png".format(K, lutnum))
    
    plt.figure()
    plt.title("Stochasticity and Effective Connectivity, LUT={0}".format(lut))
    plt.xlabel("Error rate")
    plt.ylabel("Deviation between stochastic and deterministic k_e")
    plt.scatter(err_lst2, np.array(ke_lst)-np.array(ke_dets), color='brown')
    plt.savefig("error_rate/K{0}_kedev_lutnum{1}.png".format(K, lutnum))
    
    plt.figure()
    plt.title("Stochasticity and Input Symmetry, LUT={0}".format(lut))
    plt.xlabel("Error rate")
    plt.ylabel("Deviation between stochastic and deterministic k_s")
    plt.scatter(err_lst2, np.array(ks_lst)-np.array(ks_dets), color='gray')
    plt.savefig("error_rate/K{0}_ksdev_lutnum{1}.png".format(K, lutnum))
    
#    plt.figure()
#    plt.hist(err_lst2, bins=50, color='blue')
#    plt.xlabel("Error rate")
#    plt.ylabel("Count")
#    plt.savefig("error_rate/K{0}_distr_err_lutnum{1}.png".format(K, lutnum))
    
    plt.figure()
    #Color and label quadrants for classes
    plt.scatter(kr_lst, ks_lst, color='g')
    plt.axvline(x=K/2, linestyle='-', color='r')
    plt.axhline(y=K/2, linestyle='-', color='r')
    plt.xlabel("k_r")
    plt.ylabel("k_s")
    plt.title("Input Redundancy vs. Input Symmetry, LUT={0}".format(lut))
    plt.tight_layout()
    plt.savefig("error_rate/K{0}_distr_kr_ks_lutnum{1}.png".format(K, lutnum))