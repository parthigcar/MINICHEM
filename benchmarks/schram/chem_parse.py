'''
File: chem_parse.py
Project: MINICHEM
File Created: Tuesday, 7th January 2020 8:49:49 am
Author: John Arul & Parth (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Last Modified: Thursday, 9th January 2020 8:23:33 am
Modified By: John Arul & Parth Patel (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Copyright: IGCAR - 2020
'''

import re
import numpy as np
import math
import re
from itertools import chain, combinations, permutations
from chem_parse import *


def chemparse(species, chem_name):
    """
    Parses the element and check whether all the element of the species are contained in the chem_name?
    If yes then true is returned else the false is returned.
    Ex,
    species= H2O
    chem_name = OH
    then it will return true (since both H and O are there in the chem_name)
    """

    # name = []
    # for i in range(len(chem_name)):
    if re.findall('[\+|\-]', chem_name) != []:
        return False
    elif re.findall('\([A|B|H|L|a|b|c|d|d\'|e|X|cr|I|III|II|s1|s2]*\)',
                    chem_name) != []:
        name = re.sub('\([A|B|H|L|a|b|c|d|d\'|e|X|cr|I|III|II|s1|s2]*\)',
                      "", chem_name)
    else:
        name = chem_name
    # print(name,'name')
    # chem_list = []

    chem_list = re.findall('([A-Z][a-z]*)', name)
    # print(chem_list)
    # print(chain.from_iterable(chem_list))
    chemset = set(chem_list)
    # print(chemset)
    # print(species)
    if chemset <= set(species):
        # print(chemset, set(species))
        # print('j')
        return True
    else:
        return False
