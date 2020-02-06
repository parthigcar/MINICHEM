'''
File: rf.py
Project: MINICHEM
File Created: Tuesday, 7th January 2020 12:42:07 pm
Author: John Arul & Parth (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Last Modified: Tuesday, 7th January 2020 12:49:16 pm
Modified By: John Arul & Parth Patel (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Copyright: IGCAR - 2020
'''

import re


def rf(y, species, input1, stoichiometric_dict, el_inventory):
    """
    Takes the output mole number array and returns the dictionary with the
    release fraction in the cover gas.
    input:
    y: output array containing mole number
    species: list of the species considered
    input1: list of the element initially considered.
    stoichiometric_dict: dictionary containing the stoichiometric information
    for the all the chemical species.
    el_inventory: dictionary containing the information about the input
    inventory specified.
    returns:
    - Prints the cover gas release fractions
    - writes the output in the iom.txt, released_mole_el.txt, released_sp.txt
    """
    condense_phase_match = '\([A|B|H|L|a|b|c|d|d\'|e|X|cr|I|III|II|s1|s2]*\)'
    sp_dict = {}
    released_mole_dict = {}

    considered_species = {}
    released_species = {}

    for i in range(y.shape[0]):
        sp_dict[species[i]] = y[i]

    for j in input1:
        sum2 = 0
        considered_species[j] = []
        released_species[j] = []
        for i in range(y.shape[0]):
            if re.findall(condense_phase_match, species[i]) == []:
                for i1 in stoichiometric_dict[species[i]]:
                    if i1[0].lower().capitalize() == j:
                        sum2 = sum2 + sp_dict[species[i]] * i1[1]
                        released_species[j].append(species[i])
            for i2 in stoichiometric_dict[species[i]]:
                if i2[0].lower().capitalize() == j:
                    considered_species[j].append(species[i])
        released_mole_dict[j] = sum2

    rf_dict = {}

    for i in released_mole_dict:
        rf_dict[i] = released_mole_dict[i]/el_inventory[i]

    # consolidated results
    file = open('iom.txt', 'w')
    # header of file
    file.write(f'{"Element":10s}{"mole_inventory":18s}{"Released moles":18s}\
            {"Release_fraction":18s}\t\t\t{"Released species":150s}\
            {"Considered_species":150s}\n')

    for i in sorted(rf_dict):
        file.write(f'{i:10s}{el_inventory[i]:15.3e}\
                {rf_dict[i]*el_inventory[i]:15.6e}{rf_dict[i]:20.6e}\t\t\t\
                {str(released_species[i]):150s}\
                {str(considered_species[i]):150s}\n')
    file.close()

    # all chemical species wise results
    file = open('released_mole_el.txt', 'w')
    for i in sorted(released_mole_dict):
        file.write(f'{i}{released_mole_dict[i]:15.6e}\n')
    file.close()
    print('Release fractions:', rf_dict)

    file = open('released_sp.txt', 'w')
    for i in range(y.shape[0]):
        file.write(f'{species[i]}{y[i]:15.6e}\n')
    file.close()
