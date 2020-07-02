# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:07:30 2020

@author: Sikander
"""

from itertools import product
import sys
sys.path.append('CANA')
from cana import boolean_node as bn
import pickle as pkl


def LUT_num(det_LUT):
    i = 0
    num = 0
    for inpts in det_LUT:
        num += det_LUT[inpts]*(2**i)
        i += 1
    return num

def canalization_table(K):
    """Save dict of canalization for all LUTs for a given K"""
    outputs = list(product((0,1), repeat = 2**K))
    lut = bn.BooleanNode()
    cana_dict = {}
    for num, out in enumerate(outputs):
        lut = lut.from_output_list(out)
        kr = lut.input_redundancy(norm=False)
        ke = lut.effective_connectivity(norm=False)
        ks = lut.input_symmetry(norm=False)
        cana_dict[num] = {'kr': kr, 'ke': ke, 'ks': ks} #num is the LUT #
    pkl.dump(cana_dict, open('cana_tables/K{0}_cana_table.pkl'.format(K), 'wb'))

if __name__=='__main__':
    Klist = [4]
    for k in Klist:
        canalization_table(k)
