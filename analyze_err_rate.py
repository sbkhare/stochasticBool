# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 08:21:20 2020

@author: Sikander
"""

import probabilisticBoolean as pb
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import pickle as pkl

N = 4 #sys.argv[1]
K = 2 #sys.argv[2]

def err_rate(lutnum, incr):
    err_lst = np.arange(0, 0.5, incr)
    rng = N*math.factorial(K)*(math.factorial(N-1)/(math.factorial(K)*math.factorial(N-1-K)))**N
    
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
    coh_dev_lst = []
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
            dev1, dev2, det_coherency, stoch_coherency, coh_dev = pb.stg_deviation(N, pbn, stg_det, stg_prob)
            dev1_lst.append(dev1)
            dev2_lst.append(dev2)
            coh_dev_lst.append(coh_dev)
            err_lst2.append(error)
    return p_max_lst, p_det_lst, kr_lst, ke_lst, ks_lst, classes, kr_dets, ke_dets, ks_dets, dev1_lst, dev2_lst, coh_dev_lst, err_lst2

def unique_cana(inp_num):
    ct = pkl.load(open('cana_tables/K{0}_cana_table.pkl'.format(inp_num), 'rb'))
    r = {}
    s = {}
    rs ={}
    for num in ct:
        if ct[num]['kr'] not in r.values():
            r[num] = ct[num]['kr']
        if ct[num]['ks'] not in s.values():
            s[num] = ct[num]['ks']
        if ct[num] not in rs.values():
            rs[num] = ct[num]
    return r, s, rs

if __name__=='__main__':
    lutnums = np.arange(2**(2**K))
    luts = ['C', 'NOR', '2', 'NC1', '4', 'NC2', 'XOR', 'NAND', 'AND', 'NXOR', 'C2', '11', 'C1', '13', 'OR', 'T']
    
    pm_dct = {}
    pd_dct = {}
    kr_dct = {}
    ke_dct = {}
    ks_dct = {}
    clas_dct = {}
    krd_dct = {}
    ked_dct = {}
    ksd_dct = {}
    dev1_dct = {}
    dev2_dct = {}
    coh_dev_dct = {}
    err_dct = {}
    
    for num in lutnums:
        pm, pd, kr, ke, ks, clas, krd, ked, ksd, dev1, dev2, cd, err = err_rate(num, incr=0.01)
        pm_dct[num] = pm
        pd_dct[num] = pd
        kr_dct[num] = kr
        ke_dct[num] = ke
        ks_dct[num] = ks
        clas_dct[num] = clas
        krd_dct[num] = krd
        ked_dct[num] = ked
        ksd_dct[num] = ksd
        dev1_dct[num] = dev1
        dev2_dct[num] = dev2
        coh_dev_dct[num]  = cd
        err_dct[num] = err
    r, s, rs = unique_cana(K)
    
    plt.figure()
    plt.title("Stochasticity and State Transition Graphs, N={0} K={1}".format(N, K))
    plt.xlabel("Error rate")
    plt.ylabel("Deviation between stochastic and deterministic STG edges")
    for num in lutnums:
        if num in [15]:
            plt.scatter(err_dct[num], dev1_dct[num], color="green") #label=luts[num],
    plt.legend()
    plt.savefig("error_rate/N{0}_K{1}_stgdev1.png".format(N, K))
    plt.figure()
    plt.title("Stochasticity and State Transition Graphs, N={0} K={1}".format(N, K))
    plt.xlabel("Error rate")
    plt.ylabel("Deviation between stochastic and deterministic STG non-edges")
    for num in lutnums:
        if num in [15]:
            plt.scatter(err_dct[num], dev2_dct[num], color="orange") #label=luts[num]
    plt.savefig("error_rate/N{0}_K{1}_stgdev2.png".format(N, K))
    
    plt.figure()
    plt.title("Stochasticity and Input Redundancy")
    plt.xlabel("Error rate")
    plt.ylabel("Deviation between stochastic and deterministic $k_{r}$")
    for num in lutnums:
        if num in rs:
            plt.scatter(err_dct[num], np.array(kr_dct[num])-np.array(krd_dct[num]), label="Det. $k_{r}$ = "+str(krd_dct[num][0]))
    plt.legend()
    plt.savefig("error_rate/K{0}_krdev.png".format(K))
    plt.figure()
    plt.title("Stochasticity and Effective Connectivity")
    plt.xlabel("Error rate")
    plt.ylabel("Deviation between stochastic and deterministic $k_{e}$")
    for num in lutnums:
        if num in rs:
            plt.scatter(err_dct[num], np.array(ke_dct[num])-np.array(ked_dct[num]), label="Det. $k_{e}$ = "+str(ked_dct[num][0]))
    plt.legend()
    plt.savefig("error_rate/K{0}_kedev.png".format(K))
    plt.figure()
    plt.title("Stochasticity and Input Symmetry")
    plt.xlabel("Error rate")
    plt.ylabel("Deviation between stochastic and deterministic $k_{s}$")
    for num in lutnums:
        if num in rs:
            plt.scatter(err_dct[num], np.array(ks_dct[num])-np.array(ksd_dct[num]), label="Det. $k_{s}$ = "+str(ksd_dct[num][0]))
    plt.legend()
    plt.savefig("error_rate/K{0}_ksdev.png".format(K))
    
    plt.figure()
    for num in lutnums:
        if num in rs:
            plt.scatter(kr_dct[num], ks_dct[num], label=" Det. $k_{r}$ = " +str(krd_dct[num][0]) + " $k_{s}$ = "+str(ksd_dct[num][0]))
#    plt.legend()
    plt.axvline(x=K/2, linestyle='-', color='r')
    plt.axhline(y=K/2, linestyle='-', color='r')
    plt.xlabel("$k_{r}$")
    plt.ylabel("$k_{s}$")
    plt.title("Input Redundancy vs. Input Symmetry")
    plt.tight_layout()
    plt.savefig("error_rate/K{0}_distr_kr_ks.png".format(K))
    
    plt.figure()
    for num in lutnums:
        plt.scatter(err_dct[num], coh_dev_dct[num], label=luts[num])
    plt.title("Stochasticity and Basin Coherence")
    plt.xlabel("Error rate")
    plt.ylabel("Average change in basin coherence")
    plt.savefig("error_rate/N{1}_K{0}_coherence.png".format(K, N))