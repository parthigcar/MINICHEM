'''
File: species_search.py
Project: MINICHEM
File Created: Tuesday, 7th January 2020 1:16:30 pm
Author: John Arul & Parth (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Last Modified: Tuesday, 7th January 2020 1:16:34 pm
Modified By: John Arul & Parth Patel (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Copyright: IGCAR - 2020
'''
import re
import itertools
from chem_parse import chemparse


def el(species):
    """
    The function takes the input species (list form), and convert
    the element of the list which can be compound/element to the element.
    INPUT: species (list)
    """
    sp_set = set()
    temp_list = []
    name = []
    for i in range(len(species)):
        if re.findall('\([L|a|cr|I|III|II]\)', species[i]) != []:
            name.append(re.sub('\([L|a|cr|I|III|II]\)', "", species[i]))
        # print(name)
        else:
            name.append(species[i])
    # print(name)
    for i in name:
        # print(re.findall('([A-Z][a-z]*)',i))
        temp_list.append(re.findall('([A-Z][a-z]*)', i))
        tb = list(itertools.chain.from_iterable(temp_list))
    for i in tb:
        sp_set.add(i)
    return list(sp_set)


def combination_search(species, grt_dict, combination_sp):
    """
    Searches the combination of the element in the grt_dict.
    For example,
    If species is H2O,
    Then the combinations in grt_dict might be: OH, H2O2, H, O2 etc.
    """
    a = el(species)
    for i in grt_dict.keys():
        spice = set()
        if chemparse(a, i):
            spice.add(i)
            # print(i)
            combination_sp[i] = grt_dict[i]
    # print(spice)
    return combination_sp
