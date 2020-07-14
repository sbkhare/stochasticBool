# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 19:04:04 2020

@author: Sikander
"""
import numpy as np
import random
from itertools import product
import networkx as nx
import matplotlib.pyplot as plt
import pickle as pkl

#K2_LUT_canalization = {0: {'kr': 2, 'ks': 2},
#                       1: {'kr': 0.75, 'ks': 1.5},
#                       2: {'kr': 0.75, 'ks': 0},
#                       3: {'kr': 1, 'ks': 0},
#                       4: {'kr': 0.75, 'ks': 0},
#                       5: {'kr': 1, 'ks': 0},
#                       6: {'kr': 0, 'ks': 1},
#                       7: {'kr': 0.75, 'ks': 1.5},
#                       8: {'kr': 0.75, 'ks': 1.5},
#                       9: {'kr': 0, 'ks': 1},
#                       10: {'kr': 1, 'ks': 0},
#                       11: {'kr': 0.75, 'ks': 0},
#                       12: {'kr': 1, 'ks': 0},
#                       13: {'kr': 0.75, 'ks': 0},
#                       14: {'kr': 0.75, 'ks': 1.5},
#                       15: {'kr': 2, 'ks': 2}}

def LUT_num(det_LUT): #number from LUT
    """Takes LUT as an argument, returns LUT #"""
    i = 0
    num = 0
    for inpts in det_LUT:
        num += det_LUT[inpts]*(2**i)
        i += 1
    return num

def num_LUT(K, num): #LUT from number
    """Takes LUT # as an argument, returns LUT"""
    lut = {}
    inputs = list(product((0,1), repeat = K))
    bin_num = bin(num)
    bin_num = bin_num[2:]
    if bin_num != len(inputs):
        bin_num = '0'*(len(inputs) - len(bin_num)) + bin_num
    for i, inp in enumerate(inputs):
        lut[inp] = int(bin_num[-(i+1)])
    return lut
    

#def canalizationMetrics(LUT): #For K=2 only
#    LUT_prob = {i:1 for i in range(2**(2**2))}
#    
#    for LUTnum in LUT_prob:
#        if LUTnum % 2 == 0: #if even
#            LUT_prob[LUTnum] *= LUT[(0, 0)]
#        else:
#            LUT_prob[LUTnum] *= (1 - LUT[(0, 0)])
#        
#        if LUTnum % 4 == 0 or LUTnum % 4 == 1: #if LUTnum in [0,1,4,5,8,9,12,13]
#            LUT_prob[LUTnum] *= LUT[(0, 1)]
#        else:
#            LUT_prob[LUTnum] *= (1 - LUT[(0, 1)])
#            
#        if 0 <= LUTnum <= 3 or 8 <= LUTnum <= 11:
#            LUT_prob[LUTnum] *= LUT[(1, 0)]
#        else:
#            LUT_prob[LUTnum] *= (1 - LUT[(1, 0)])
#        
#        if LUTnum <= 7:
#            LUT_prob[LUTnum] *= LUT[(1, 1)]
#        else:
#            LUT_prob[LUTnum] *= (1-LUT[(1, 1)])
#    
##    for LUTnum in LUT_prob:
##        bin_num = bin(LUTnum)
##        bin_num = bin_num[2:]
##        for idx, inp in enumerate(LUT):
##            if bin_num[-(idx+1)] == '0':
##                LUT_prob[LUTnum] *= LUT[inp]
##            elif bin_num[-(idx+1)] == '1':
##                LUT_prob[LUTnum] *= (1 - LUT[inp])
#            
#    kr = 0
#    ks = 0
#    for LUTnum in LUT_prob:
#        kr += LUT_prob[LUTnum] * K2_LUT_canalization[LUTnum]['kr']
#        ks += LUT_prob[LUTnum] * K2_LUT_canalization[LUTnum]['ks']
#    ke = 2 - kr
#    clas = 'str'
#    if kr >= 1 and ks >= 1:
#        clas = 'A'
#    elif kr >= 1 and ks < 1:
#        clas = 'B'
#    elif kr < 1 and ks >= 1:
#        clas = 'C'
#    elif kr < 1 and ks < 1:
#        clas = 'D'
#    return kr, ke, ks, clas, LUT_prob

def sortvals(dct, order=False):
    """Sort values of a dict in descending order"""
    nd = dict()
    for k in sorted(dct, key=dct.get, reverse=order):
        nd[k] = dct[k]
    return nd

def saveGML(G, name):
    mapping = {node: str(node) for node in G.node}
    G = nx.relabel_nodes(G, mapping)
    nx.write_gml(G, name)

def absoluteThreshold(H, thresh, show=False): #Removes edges in STG with probability less than threshold
    """Removes edges in State Transiton Graph with probability less than threshold, H = state transition graph, thresh = probability threshold"""
    G = H.copy()
    remove = []
    for trans in G.edges():
        if G[trans[0]][trans[1]]['probability'] < thresh:
            remove.append(trans)
    G.remove_edges_from(remove)
    if show:
        plt.figure()
        nx.draw_circular(G)
        nx.draw_networkx_labels(G, pos=nx.circular_layout(G))
        nx.draw_networkx_edge_labels(G, pos=nx.circular_layout(G))
    return G

def cumulativeThreshold(H, thresh, show=False):
    """Keeps edges the minimum number of edges in State Transiton Graph needed 
    to have an outgoing probability from a state be higher than the threshold, 
    H = state transition graph, thresh = probability threshold"""
    G = H.copy()
    remove = []
    for node in G.nodes():
        cum_p = 0
        edge_prob_dict = {(n1, n2): stg_prob[n1][n2]['probability'] for (n1, n2) in stg_prob.edges(node)}
        epd = sortvals(edge_prob_dict, order=True)
        for edge in epd:
            if cum_p < thresh:
                cum_p += epd[edge]
            else:
                remove.append(edge)
    G.remove_edges_from(remove)
    plt.figure()
    nx.draw_circular(G)
    nx.draw_networkx_labels(G, pos=nx.circular_layout(G))
#    nx.draw_networkx_edge_labels(G, pos=nx.circular_layout(G))
    return G

def stg_deviation(N, sbn, det_stg, stoch_stg):
    """Calculates deviation between deterministic and stochastic state transition graphs"""
    dev1 = 0 #Deterministic-transitions
    dev2 = 0 #Error-transitions
    i = 0
    for s_edge in stoch_stg.edges.data():
        i+=1
        if det_stg.has_edge(s_edge[0], s_edge[1]):
            dev1 += 1 - s_edge[2]['probability']
        else:
            dev2 += s_edge[2]['probability']
    dev1 = dev1/(2**N) #D-TRANSITIONS
    dev2 = dev2/((2**N)**2 - 2**N) #E-TRANSITIONS
    print("State transition graph deviation: ", dev1, ", ", dev2)
    #ATTRACTOR ANALYSIS
    coherency_det = {}
    coherency_stoch = {}
    for attr in sbn.attractors_det:
        #JUST DO OBSERVED COHERENCY
        basin_coh_det = 0
        basin_coh_prob = 0
        basin_states = set()
        #Find set ofbasin states
        for nd in det_stg.nodes():
            if nx.has_path(det_stg, source=nd, target=attr[0]):
                basin_states.add(nd)
        #Calculate coherency of each of the basin states
        for bs in basin_states:
            state_coh_det = 0
            state_coh_prob = 0
            for i in range(sbn.N):
                #Perturb one of the node's state
                perturb_state = list(bs)
                if perturb_state[i] == 1:
                    perturb_state[i] = 0
                else:
                    perturb_state[i] = 1
                perturb_state = tuple(perturb_state)
                if perturb_state in basin_states:
                    state_coh_det += 1
                    for bs1 in basin_states:
                        state_coh_prob += stoch_stg[perturb_state][bs1]['probability']
                state_coh_det /= sbn.N
                state_coh_prob /= sbn.N
                basin_coh_det += state_coh_det
                basin_coh_prob += state_coh_prob
        basin_coh_det /= len(basin_states) + 1
        basin_coh_prob /= len(basin_states) + 1
        if len(attr) == 1:
            coherency_det[attr[0]] = basin_coh_det
            coherency_stoch[attr[0]] = basin_coh_prob
        elif len(attr) > 1:
            coherency_det[tuple(attr)] = basin_coh_det
            coherency_stoch[tuple(attr)] = basin_coh_prob
    coherency_dev = 0
    for attr in coherency_det:
        coherency_dev +=  coherency_stoch[attr] - coherency_det[attr]
    coherency_dev /= len(coherency_det)
    return dev1, dev2, coherency_det, coherency_stoch, coherency_dev #attr_leave_prob, avg_alp
    
#Also need function for keeping edges that account for x% of probability 

class PBNetwork():
    def __init__(self, N, K):
        """NK network in which all nodes have the same LUT"""
        self.N = N
        self.K = K
        self.state = np.zeros(N, dtype=int)
        self.adjmat = np.zeros((N, N))
        self.LUT = {} #values are p0
        self.LUT_det = {}
        self.LUT_alt = {}
        self.ct = pkl.load(open('cana_tables/K{0}_cana_table.pkl'.format(K), 'rb'))
        self.attractors_det = []
        
    def randomizeInitState(self):
        """Randomize initial state"""
        self.state = np.random.randint(2, size=self.N)
        print("0: ", self.state)
        
    def createNetwork(self):
        """Creates network (no self loops)"""
        rownum = 0
        for row in self.adjmat:
            ind = [i for i in range(self.N)]
            ind.remove(rownum) #No self loops
            sample = random.sample(ind, self.K)
            for i in range(self.K):
                row[sample[i]] = 1
            rownum += 1
        
    def randomizeLUT(self):
        """Chooses a random LUT for the given K. Probaiblities are chosen from a uniform distribution from 0 to 1"""
        inputs = list(product((0,1), repeat = self.K))
        for i in inputs:
            self.LUT[i] = np.random.rand()
        print("Randomized LUT:")
        print(self.LUT)
        
    def biasedLUT(self, det=True):
        """Chooses an LUT that is either biased towards being more determinsitic or stochastic"""
        inputs = list(product((0,1), repeat = self.K))
        if det: #biased towards deterministic, CHANGE SO THAT ONLY AVERAGE P-MAX > 0.75 IS ACCEPTED
            for i in inputs:
                self.LUT[i] = np.random.beta(0.5, 0.5) #Probabilities drawn from beta distribution which peaks at 0 and 1
            print("Deterministic biased LUT:")
            print(self.LUT)
        else: #biased towards probabilistic
            for i in inputs:
                p0 = np.random.normal(0.5, 0.18) #Probabilities drawn from a normal distribution with mu=0.5
                while p0 > 1 or p0 < 0:
                    p0 = np.random.normal() 
                self.LUT[i] = p0
            print("Probabilistic biased LUT:")
            print(self.LUT)
            
    def biasedRangeLUT(self, mn, mx):
        """Chooses an LUT with a pmax within the range [mn, mx]"""
        #pmax range is biased drectly
        inputs = list(product((0,1), repeat = self.K))
        for i in inputs:
            if np.random.rand() < 0.5:
                self.LUT[i] = np.random.uniform(mn, mx)
            else:
                self.LUT[i] = np.random.uniform(1-mx, 1-mn)
        print("Biased range LUT: {0} < p < {1}".format(mn, mx))
        print(self.LUT)
    #Create functions for biasing towards and, or, xor
            
    def errLUT(self, lutnum, err):
        """Chooses an LUT where each entry has an an erro rate, lutnum = deterministic LUT #, err = error rate"""
        lut = num_LUT(self.K, lutnum)
        for inp in lut:
            if lut[inp] == 0:
                self.LUT[inp] = 1 - err
            else:
                self.LUT[inp] = err
        print("Most likely deterministic LUT: CANA #{0}".format(lutnum))
        print(self.LUT)
    
    
    def deterministicLUT(self):
        """find most likely LUT and its probability"""
        p_det = 1
        p = 0
        for inpt in self.LUT:
            if self.LUT[inpt] < 0.5: #if p0 <0.5 then it should be 1 in the most likely deterministic case
                self.LUT_det[inpt] = 1
                p += (1 - self.LUT[inpt])
                p_det *= (1 - self.LUT[inpt])
            else:
                self.LUT_det[inpt] = 0
                p += self.LUT[inpt]
                p_det *= self.LUT[inpt]
        p = p/(2**self.K)
        num = LUT_num(self.LUT_det)
        print("Most likely deterministic LUT: CANA #{0}".format(num))
        print(self.LUT_det)
        print("Average p_max = {0}".format(p))
        print("Probability of most likely deterministic LUT:", p_det)
        kr = self.ct[num]['kr']
        ke = self.ct[num]['ke']
        ks = self.ct[num]['ks']
        return p, p_det, kr, ke, ks
        
        #return Probability

    def alternativeLUT(self):
        """find alternative LUT and its probability"""
        p = 0
        for inpt in self.LUT:
            if self.LUT[inpt] < 0.5: #if p0 > 0.5 then it should be 0 in the most likely deterministic case
                self.LUT_alt[inpt] = 0
                p += self.LUT[inpt]
            else:
                self.LUT_alt[inpt] = 1
                p += (1 - self.LUT[inpt])
        p = p/(2**self.K)
        print("Alternative LUT: average p_min = {0}".format(p))
        print(self.LUT_alt)
        return p
        #return probability
    
    def canalization(self, pmax, show=False, save=False):
        """Calculates stochastic canlization"""
        LUT_prob = {i:1 for i in range(2**(2**self.K))}
        for LUTnum in LUT_prob:
            bin_num = bin(LUTnum)
            bin_num = bin_num[2:]
            if bin_num != len(self.LUT):
                bin_num = '0'*(len(self.LUT) - len(bin_num)) + bin_num
            for idx, inp in enumerate(self.LUT):
                if bin_num[-(idx+1)] == '0':
                    LUT_prob[LUTnum] *= self.LUT[inp]
                elif bin_num[-(idx+1)] == '1':
                    LUT_prob[LUTnum] *= (1 - self.LUT[inp])
        kr = 0
        ke = 0
        ks = 0
        for LUTnum in LUT_prob:
            kr += LUT_prob[LUTnum] * self.ct[LUTnum]['kr']
            ke += LUT_prob[LUTnum] * self.ct[LUTnum]['ke']
            ks += LUT_prob[LUTnum] * self.ct[LUTnum]['ks']
        clas = 'str'
        if kr >= 1 and ks >= 1:
            clas = 'A'
        elif kr >= 1 and ks < 1:
            clas = 'B'
        elif kr < 1 and ks >= 1:
            clas = 'C'
        elif kr < 1 and ks < 1:
            clas = 'D'
        if show:
            plt.figure()
            plt.bar(LUT_prob.keys(), LUT_prob.values(), color='g')
            plt.title("Distribution of Look-up Tables, p_max = {0}".format(round(pmax, 4)))
            plt.ylabel("Probability")
            plt.xlabel("LUT #")
            plt.ylim(0, max(LUT_prob.values()))
            plt.show()
        if save:
            plt.savefig("NK/N{0}_K{1}_pmax{2}_lutdistr.png".format(self.N, self.K, pmax))
        return kr, ke, ks, clas, LUT_prob
        
    def nextState(self):
        """Iterates to the next state"""
        newstate = np.zeros(self.N, dtype=int)
        i = 0
        for adj in self.adjmat:
            input_ind = np.where(adj == 1)
            inputs = [self.state[ind] for ind in input_ind[0]]
            if np.random.rand() < self.LUT[tuple(inputs)]:
                newstate[i] = 0
            else:
                newstate[i] = 1
            i += 1
        return newstate
        
    def run(self, numsteps):
        """Random walk with specified number of steps, returns STG with number 
        of times an edge was traversed as the edge weight and number of times 
        A STATE WAS VISITED AS THE NODE ATTRIBUTE"""
        state_trans_G = nx.DiGraph()
        state_trans_G.node[tuple(self.state)]['count'] = 1
        for i in range(1, numsteps):
            newstate = self.nextState()
            if state_trans_G.has_edge(tuple(self.state), tuple(newstate)):
                state_trans_G[tuple(self.state)][tuple(newstate)]['count'] += 1
                state_trans_G.node[tuple(newstate)]['count'] = 1
            else:
                state_trans_G.add_edge(tuple(self.state), tuple(newstate), count=1)
                state_trans_G.node[tuple(newstate)]['count'] = 1
            self.state = newstate
            print("{0}: ".format(i), self.state)
        nx.draw(state_trans_G, pos=nx.shell_layout(state_trans_G))
        nx.draw_networkx_labels(state_trans_G, pos=nx.shell_layout(state_trans_G))
        nx.draw_networkx_edge_labels(state_trans_G, pos=nx.shell_layout(state_trans_G))
        return state_trans_G
        
    def stateTransitionGraph(self, pmax, show=False, save=False):
        """Returns the stochastic state transition graph with transition probability as edge weights"""
        G = nx.DiGraph()
        G.add_nodes_from(list(product((0,1), repeat = self.N)))
        plist = []
        for state1 in G.nodes:
            p0 = []
            for adj in self.adjmat:
                input_ind = np.where(adj == 1)
                inputs = [state1[ind] for ind in input_ind[0]]
                p0.append(self.LUT[tuple(inputs)])
            for state2 in G.nodes:
                p_edge = 1
                i = 0
                for s1, s2 in zip(state1, state2):
                    if s2 == 0:
                        p_edge *= p0[i]
                    else:
                        p_edge *= (1 - p0[i])
                    i += 1
                p_edge = round(p_edge, 4)
                plist.append(p_edge)
                G.add_edge(state1, state2, probability=p_edge)
        if show:
            plt.figure()
            nx.draw_circular(G)
            nx.draw_networkx_labels(G, pos=nx.circular_layout(G))
#            nx.draw_networkx_edge_labels(G, pos=nx.circular_layout(G))
            print("Mean edge probability: ", np.mean(plist))
            plt.figure()
            plt.hist(plist, bins=2**self.N, color='orange') #, log=True
            plt.title("Edge probability distribution, p_max = {0}".format(round(pmax, 4)))
            plt.xlabel("Edge probability")
            plt.ylabel("Count")
            plt.show()
        if save:
            plt.savefig("NK/N{0}_K{1}_pmax{2}_edgprob.png".format(self.N, self.K, pmax))
        return G, np.mean(plist)+np.std(plist)
    
    def detSTG(self, show=False):
        """Returns the determistic state transition graph"""
        G = nx.DiGraph()
        G.add_nodes_from(list(product((0,1), repeat = self.N)))
        for state in G.nodes:
            newstate = []
            for adj in self.adjmat:
                input_ind = np.where(adj == 1)
                inputs = [state[ind] for ind in input_ind[0]]
                newstate.append(self.LUT_det[tuple(inputs)])
            G.add_edge(state, tuple(newstate))
#            if state == tuple(newstate):
#                self.attractors_det.append([state])
        self.attractors_det = list(nx.simple_cycles(G))
        if show:
            plt.figure()
            nx.draw_circular(G)
            nx.draw_networkx_labels(G, pos=nx.circular_layout(G))
        return G     
                
if __name__=='__main__':
    Ni = 3
    Ki = 2
    pbn = PBNetwork(N=Ni, K=Ki) #Number of possible combinations of deterministic LUTs and adjacency matrixes: (2**(2**K))*(math.factorial(N-1)/math.factorial(N-1-K))**N
    pbn.createNetwork()
##    pbn.randomizeLUT()
#    pbn.biasedLUT(det=True)
#    pbn.biasedRangeLUT(mn=0.5, mx=0.55)
    pbn.errLUT(10, 0.25)
    pmax, pdet, kr, ke, ks = pbn.deterministicLUT()
    kr, ke, ks, clas, LUT_distr = pbn.canalization(pmax=pmax, show=True, save=False)
    print("Class ", clas)
    print("k_r = ", kr)
    print("k_e = ", ke)
    print("k_s = ", ks)
##    pbn.randomizeInitState()
    stg_prob, p_thresh = pbn.stateTransitionGraph(pmax=pmax, show=True, save=False)
##    print("STG edges and probabilities")
##    for edge in stg_prob.edges.data():
##        print(edge)
#    
##    pbn.alternativeLUT()
    stg_det = pbn.detSTG(show=True)
##    print("STG LUT edges:")
##    print(stg_det.edges())
    dev1, dev2, det_coherency, stoch_coherency, coh_dev = stg_deviation(Ni, pbn, stg_det, stg_prob)
    stg_thresh = absoluteThreshold(stg_prob, thresh=(0.5)**Ni, show=True) #thresh=0.5 for deterministic biased and = 0.15 for probabilistic biased
#    stg_thresh = cumulativeThreshold(stg_prob, thresh=0.6, show=True)
#    pmax = round(pmax, 4)
#    pmax = str(pmax)
#    p_thresh = round(p_thresh, 4)
#    p_thresh = str(p_thresh)
#    nx.write_graphml(stg_det, "NK/N{0}_K{1}_pmax{2}_stgdet.graphml".format(Ni, Ki, pmax[2:]))
#    nx.write_graphml(stg_thresh, "NK/N{0}_K{1}_pmax{2}_stgthresh__thr{3}.graphml".format(Ni, Ki, pmax[2:], p_thresh[2:]))
#    pkl.dump(pbn.LUT, open("NK/N{0}_K{1}_pmax{2}_lut.pkl".format(Ni, Ki, pmax[2:], p_thresh[2:]), 'wb'))

    
    