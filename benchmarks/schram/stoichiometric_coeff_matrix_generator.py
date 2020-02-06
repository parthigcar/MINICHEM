'''
File: stoichiometric_coeff_matrix_generator.py
Project: MINICHEM
File Created: Tuesday, 7th January 2020 11:12:49 am
Author: John Arul & Parth (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Last Modified: Thursday, 9th January 2020 8:26:47 am
Modified By: John Arul & Parth Patel (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Copyright: IGCAR - 2020
'''

import numpy as np


def stoi(species, input1, stoichiometric_dict):
    """
    Description:
    Makes one row of the stoichiometry coefficient.

    Parameters:
    species:List of species
    input1: list of the elements considered.
    stoichiometric_dic: dictionary of the species with the stoichiometric
                        information

    Returns:
    list of the row of the stoichiometric matrix
    """

    list1 = np.zeros(len(input1))
    for i in range(len(input1)):
        for j in stoichiometric_dict[species]:
            if input1[i] == j[0].lower().capitalize():
                # print(input1[i], j[0])
                list1[i] = j[1]
    return list1


def make_ac(input1, b, considered_sp_c, stoichiometric_dict):
    """
    Description:
    Makes stoichiometry matrix for the condensed species.

    Parameters:
    input1: list of the elementts
    b: element inventory value
    considered_sp_c: set of considered species in the condensed phase

    Returns:
    a_c: condensed stoichiometry matrix
    """
    a_c = np.zeros([len(b), len(considered_sp_c)])
    k2 = 0
    for j in considered_sp_c:
        a_c[:, k2] = stoi(j, input1, stoichiometric_dict)
        k2 = k2 + 1
    return a_c


def test_for_dependence(a_c, inds, input1, b, stoichiometric_dict,
                        sp_c, dict_of_all_sp_grt, a, pis, initial_sp_c,
                        total_sp_c):
    """
    Checks the a_c matrix for the dependent row, and returns the list of the
    dependent rows.
    """
    list_pot_removed = []
    list_of_dependent_sp = []
    k34 = 0
    for i in range(a_c.transpose().shape[0]):
        if i not in inds:
            # print(i)
            # dependent_sp_c = sp_c[i]
            list_of_dependent_sp.append(sp_c[k34])
            # total_sp_c.append(sp_c[i])
            print('Species', sp_c[k34], 'Removed from sp_c list')
            del sp_c[k34]
            k34 = k34-1
        k34 = k34 + 1

    # a_c = make_ac(input1, b, sp_c, stoichiometric_dict)
    # print('The dependent sp c is:', dependent_sp_c, sp_c)
    # list_of_dependent_sp.append(dependent_sp_c)
    temp3 = sp_c.copy()
    s1 = set()
    for dependent_sp_c in list_of_dependent_sp:
        for i in stoichiometric_dict[dependent_sp_c]:
            s1.add(i[0])
    s1.discard('*')
    k14 = 0
    for i in range(len(temp3)):
        s2 = set()
        for i1 in stoichiometric_dict[temp3[i]]:
            s2.add(i1[0])
        s2.discard('*')
        s2.discard('O')
        if s2 <= s1:
            # print(s1, s2)
            list_of_dependent_sp.append(sp_c[k14])
            # total_sp_c.append(sp_c[k14])
            print('Species removed:', sp_c[k14])
            del sp_c[k14]
            k14 = k14 - 1
        k14 = k14 + 1
    print('List of dependent species:', list_of_dependent_sp)
    for i in range(len(list_of_dependent_sp)):
        try:
            del sp_c[sp_c.index(list_of_dependent_sp[i])]
            total_sp_c.append(list_of_dependent_sp[i])
            print('Deleted sp_c', list_of_dependent_sp[i])
            print('Added total_sp_c', list_of_dependent_sp[i])
        except:
            total_sp_c.append(list_of_dependent_sp[i])
            print('Added total_sp_c', list_of_dependent_sp[i])
    # a_c = make_ac(input1, b, sp_c, stoichiometric_dict)

        potent = (dict_of_all_sp_grt[list_of_dependent_sp[i]] -
                  sum(pis *
                      a[:, initial_sp_c.index(list_of_dependent_sp[i])]))
        print(potent, list_of_dependent_sp[i])
        # if potent < 0:
        list_pot_removed.append(potent)

    if len(list_of_dependent_sp) > 1:
        print('Added in sp_c',
              list_of_dependent_sp[
                  list_pot_removed.index(min(list_pot_removed))])
        sp_c.append(list_of_dependent_sp[
            list_pot_removed.index(min(list_pot_removed))])
        try:
            del total_sp_c[total_sp_c.index(
                           list_of_dependent_sp[
                               list_pot_removed.index(min(list_pot_removed))])]
            print('deleted from total_sp_c',
                  total_sp_c[total_sp_c.index(
                      list_of_dependent_sp[
                          list_pot_removed.index(min(list_pot_removed))])])
        except:
            # print('Passed exception in del total_sp_c')
            pass
        a_c = make_ac(input1, b, sp_c, stoichiometric_dict)
    elif len(list_of_dependent_sp) == 1:
        # total_sp_c.append(list_of_dependent_sp[0])
        sp_c.append(
            list_of_dependent_sp[list_pot_removed.index(min(list_pot_removed))])
        # print('HI')
        try:
            del total_sp_c[total_sp_c.index(
                list_of_dependent_sp[
                    list_pot_removed.index(min(list_pot_removed))])]
            print('Deleted from total_sp_c', total_sp_c[total_sp_c.index(
                list_of_dependent_sp[
                    list_pot_removed.index(min(list_pot_removed))])])
        except:
            # print('Passed exception in del total_sp_c')
            pass
        a_c = make_ac(input1, b, sp_c, stoichiometric_dict)
    for i in list_of_dependent_sp:
        if i != list_of_dependent_sp[
                list_pot_removed.index(min(list_pot_removed))] and i\
                not in total_sp_c:
            total_sp_c.append(i)
            print(i, 'Added back in total_sp_c')

    return a_c, sp_c, total_sp_c, list_of_dependent_sp
