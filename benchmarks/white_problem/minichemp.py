'''
File: minichemp.py
Project: MINICHEM (Minimisation of CHEMical potential)
File Created: Tuesday, 9th April 2019 12:44:07 pm
Author: John Arul & Parth (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Last Modified: Thursday, 18th April 2019 1:32:54 pm
Modified By: John Arul & Parth Patel (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Thermochemical minimization using SLSQP and gradient descent
Copyright: IGCAR - 2019
To run this program following files are required:
- chem_parse.py
- data_process.py
- mini.py
- plot_hv.py
- pythormoread.py
- rf.py
- species_search.py
- stoichiometric_coeff_matrix_generator.py
Input files:
- thermo_chemical_database.txt
- thermo_python.in
Output files:
- iom.txt: Gives the released mole, RF, and details about considered species
            and released species.
- released_sp.txt: Gives the released species mole
- released_mole_el.txt: el wise result of the released mole.
'''

import re
import rf
import sympy
import warnings
import numpy as np
from scipy import optimize
from itertools import chain

# imports from other modules
import mini
import plot_hv
import data_process
import pythermoread
import species_search
import stoichiometric_coeff_matrix_generator

# =============================================================================


# =============================================================================
# White case
input_sp = ['2H', '1N', '1O']
# ALMR Case
# input_sp = ['3.10e3He', '1.24e2Kr', '1.46e3Xe', '0.87e2O', '3.10H',
#             '1.22e2I', '0.66e1Br', '3.02e7Na', '1.24e3Cs', '1.04e2Rb',
#             '0.43e3Ba', '2.35e2Sr', '1.36e3Mo', '1.22e3Zr', '0.59e5U',
#             '2.12e1Sb']

# input compounds converted to elements
input1, b = data_process.ip_sp_to_el(input_sp)
el_inventory = {}
k23 = 0
for i in input1:
    el_inventory[i] = b[k23]
    k23 = k23 + 1

v = 140    # m3    100 for PSFR and 140 for ALMR
R = 8.31445984848484848484   # gas constant
temperature = 3500  # 883 + 273    # 90 + 273   # K
pressure = 51.7106797    # bar
method = 'MINI'   # Valide option: MINI or SLSQP
switch = 'TP'     # Valid option: TP or TV
trace = 1e-25
combination_sp = {}
grt_dict = {}
stoichiometric_dict = {}
thermo_dict, stoichiometric_dict = pythermoread.thermoread()
dict_of_all_sp_grt = \
    pythermoread.calculate_grt(grt_dict, temperature, thermo_dict)


# =============================================================================
# only GRT portion
# list_sp_according_schram = ['Am', 'Ba', 'Ba2', 'BaI', 'BaI2', 'Br', 'Br2',
#                             'Ce', 'CeO', 'CeO2', 'Cm', 'Cs', 'Cs2', 'Cs2I2',
#                             'CsI', 'Eu', 'H2', 'H2O', 'He', 'I', 'I2', 'Kr',
#                             'La', 'Mo', 'Na', 'Na2', 'Na2Br2', 'Na2I2',
#                             'NaBr', 'NaI', 'Np', 'O2', 'PuO2', 'Rb', 'Sr',
#                             'Sr2', 'SrI', 'SrI2', 'U', 'UO2', 'Xe', 'Zr',
#                             'Sb', 'Sb2', 'Sb4', 'Te', 'Te2', 'ZrI',
#                             'ZrI2', 'ZrI3', 'ZrI4', 'Ru', 'ZrBr4', 'U(L)',
#                             'Na(L)', 'Am(L)', 'UO2(cr)', 'Eu(L)', 'Sb(L)',
#                             'Np(L)', 'Ce(L)', 'La(c)', 'Pu(L)', 'Mo(L)',
#                             'Cm(b)', 'Ru(cr)', 'Zr(L)', 'Cs2Te(L)',
#                             'ZrTe2(cr)', 'RuTe(cr)', 'CsI(L)', 'Na(L)',
#                             'Ba(L)', 'CsBr(L)', 'CsI(L)', 'Ba(L)',
#                             'NaI(L)', 'Rb(L)', 'BaI2(cr)', 'SrI2(L)',
#                             'Ce2O3(cr)', 'TeI4(cr)',  'BaI2(L)', 'CsBr(L)',
#                             'NaI(L)', 'Rb(L)', 'Sr(L)', 'Cs(L)'
#                             'Cs(L)', 'Sr(L)', 'SrI2(L)', 'RuTe2(cr)',
#                             'PuO2(L)', 'ZrBr4(cr)', 'La(L)', 'Zr(b)',
#                             'CeO2(L)', 'Ru(L)', 'PuO2(cr)', 'ZrI4(cr)']

# finding combination
grt_dict1 = species_search.combination_search(
    input1, dict_of_all_sp_grt, combination_sp)
# grt_dict1 = pythermoread.only_grt(grt_dict, list_sp_according_schram)

# =============================================================================
grt_dict = {}
for i in grt_dict1:
    if grt_dict1[i] > 1000:
        pass
    else:
        grt_dict[i] = grt_dict1[i]

# Externally providing chemical potential
grt_dict = {}
grt_dict['H'] = -10.021
grt_dict['H2'] = -21.096
grt_dict['H2O'] = -37.986
grt_dict['N'] = -9.846 
grt_dict['N2'] =-28.653
grt_dict['NH'] = -18.918
grt_dict['NO'] = -28.032
grt_dict['O'] = -14.640
grt_dict['O2'] = -30.594
grt_dict['OH'] = -26.111

species = []

for i in grt_dict:
    species.append(i)


condense_phase_match = '\([A|B|H|L|a|b|c|d|d\'|e|X|cr|I|III|II|s1|s2]*\)'


if method == 'MINI':
    sp_g = []
    total_sp_c = []

    for i in range(len(species)):
        if re.findall(condense_phase_match, species[i]) == []:
            sp_g.append(species[i])
        else:
            total_sp_c.append(species[i])

    # Ag matrix generation, this matrix size is not changing
    # hence this matrix is out side of the mini module and fed
    # as input to mini module.
    a_g = np.zeros([len(b), len(sp_g)])
    k1 = 0
    for j in sp_g:
        a_g[:, k1] = stoichiometric_coeff_matrix_generator.stoi(
            j, input1, stoichiometric_dict)
        k1 = k1 + 1
    # stoichiometric matrix for solids
    a = np.zeros([len(b), len(total_sp_c)])
    k1 = 0
    for j in total_sp_c:
        a[:, k1] = stoichiometric_coeff_matrix_generator.stoi(
            j, input1, stoichiometric_dict)
        k1 = k1 + 1
    initial_sp_c = total_sp_c.copy()

    # INSERT = ['Cs(L)', 'Na(L)', 'UO2(cr)', 'PuO2(cr)', 'Ru(cr)']
    INSERT = []

    y, sp_g, sp_c = mini.mini_solver(input1, b, sp_g, INSERT, total_sp_c,
                                     a, a_g, trace, dict_of_all_sp_grt,
                                     initial_sp_c, grt_dict,
                                     stoichiometric_dict,
                                     switch, temperature, v=v,
                                     pressure=pressure)

    # =========================================================================
    #                          Release fraction module
    # =========================================================================
    species_updated = sp_g + sp_c
    rf.rf(y, species_updated, input1, stoichiometric_dict, el_inventory)

elif method == 'SLSQP':
    a = np.zeros([len(input1), len(species)])
    x0 = np.ones(len(species)) * 0.1

    k1 = 0
    for j in species:
        a[:, k1] = stoichiometric_coeff_matrix_generator.stoi(
            j, input1, stoichiometric_dict)
        k1 = k1 + 1
    no_it = 20
    opt1 = {'eps': 1e-3, 'maxiter': 20000,
            'ftol': 1e-6, 'iprint': 3, 'disp': False}
    opt2 = {'eps': 1e-3, 'maxiter': 2500,
            'ftol': 1e-6, 'iprint': 3, 'disp': False}

    def constr(x):
        return a.dot(x)-b
    scale = 1
    for i in range(no_it):
        # print(list(zip(species, x0)))
        print('Iteration number:', i)
        if i > 0:
            # print(res.message)
            k1 = 0
            species = np.array(species)
            x0 = res.x
            bnds = [(0.0, None) for i in range(x0.shape[0])]
            if switch == 'TV':
                res = optimize.minimize(mini.min_fun_helmholtz, x0,
                                        args=(species, grt_dict,
                                              temperature, v),
                                        method='SLSQP',
                                        constraints=[{'type': 'eq', 'fun':
                                                      constr}],
                                        bounds=bnds, options=opt2)
            elif switch == 'TP':
                res = optimize.minimize(mini.gibbs_calculate, x0,
                                        args=(species, grt_dict,
                                              temperature, pressure),
                                        method='SLSQP',
                                        constraints=[{'type': 'eq', 'fun':
                                                      constr}],
                                        bounds=bnds, options=opt2)
            if res.success:
                print('success')
                break
        else:
            bnds = [(0.0, None) for i in range(x0.shape[0])]
            if switch == 'TV':
                res = optimize.minimize(mini.min_fun_helmholtz, x0,
                                        args=(species, grt_dict,
                                              temperature, v),
                                        method='SLSQP',
                                        constraints=[{'type': 'eq', 'fun':
                                                      constr}],
                                        bounds=bnds, options=opt2)
            elif switch == 'TP':
                res = optimize.minimize(mini.gibbs_calculate, x0,
                                        args=(species, grt_dict,
                                              temperature, pressure),
                                        method='SLSQP',
                                        constraints=[{'type': 'eq', 'fun':
                                                      constr}],
                                        bounds=bnds, options=opt2)
            if res.success:
                print('success')
                break

    # print(list(zip(species, res.x)))
    abc = res.success
    y = res.x
    rf.rf(y, species, input1, stoichiometric_dict, el_inventory)

else:
    print('Please enter valid method')
include_el = ['N', 'O', 'H']
# (allowed formats:
# ['html', 'json', 'auto', 'png', 'widgets', 'scrubber', 'auto', None])
plot_hv.plot_hv(input1, stoichiometric_dict, include_el, 'white_problem.html', 
                Min=0, Max=1e3)
print('Program execution completed')
