'''
File: data_process.py
Project: MINICHEM
File Created: Tuesday, 7th January 2020 8:49:50 am
Author: John Arul & Parth (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Last Modified: Thursday, 9th January 2020 8:23:45 am
Modified By: John Arul & Parth Patel (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Copyright: IGCAR - 2020
'''

import numpy as np
import itertools
from itertools import chain
import re
import scipy
import warnings
from scipy.linalg import solve


def grt(g_hts, hf, T):
    """
    This function is written for comparing the grt value
    with thermo read value
    Function takes g_hts values and hf values and calculate value for g/rt
    """
    R1 = 8.31445984848484848484    # J/K/mol
    return print(-g_hts/R1 + hf/R1/T * 1000)


def ip_sp_to_el(input_sp):
    '''
    converts the input species into element. And stores the input moles of the
    particalular species.
    Input:
    input_sp: takes list of the species
    Output:
    input1: converted list of element
    b: no of moles for that element
    '''
    element = []
    el_dict = {}
    for i in input_sp:
        a = re.findall('(\d*\.?\d*[E|e]?[+|-]?\d+)*([A-Z][a-z]?)(\d)*', i)
        element.append(a)
    el_copy = list(chain.from_iterable(element))
    b = {}
    s = set()
    for i in el_copy:
        s.add(i[1])
    for i in list(s):
        b[i] = 0

    for k in element:
        if len(k) >= 2:
            if k[0][0] != '':
                multiplier = float(k[0][0])
            else:
                multiplier = 1
            for l in k:
                if l[2] == '':
                    b[l[1]] = float(l[2].replace('', '1')) * multiplier + \
                              b[l[1]]
                else:
                    b[l[1]] = float(l[2]) * multiplier + b[l[1]]
        else:
            k = k[0]
            if k[0] != '' and k[2] == '':
                b[k[1]] = float(k[2].replace('', '1'))*float(k[0]) + b[k[1]]
            elif k[2] != '' and k[0] == '':
                b[k[1]] = float(k[2]) * float(k[0].replace('', '1')) + b[k[1]]
            elif k[2] == '' and k[0] == '':
                b[k[1]] = float(k[0].replace('', '1')) * \
                          float(k[2].replace('', '1')) + b[k[1]]
            else:
                b[k[1]] = float(k[2]) * float(k[0]) + b[k[1]]
    return list(b.keys()), list(b.values())
