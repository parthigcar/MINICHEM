'''
File: mini.py
Project: MINICHEM
File Created: Tuesday, 7th January 2020 10:46:58 am
Author: John Arul & Parth (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Last Modified: Tuesday, 7th January 2020 10:54:36 am
Modified By: John Arul & Parth Patel (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Copyright: IGCAR - 2020
'''

import re
import numpy as np
import warnings
import sympy
import scipy
from scipy.linalg import solve
import stoichiometric_coeff_matrix_generator

"""
This file contains mini solver (based on the quadratic gradient descent method)
and other helping functions, which are used in the solver.
"""


def mini_solver(input1, b, sp_g, INSERT, total_sp_c, a, a_g, trace,
                dict_of_all_sp_grt, initial_sp_c, grt_dict,
                stoichiometric_dict, switch,
                temperature, v=0, pressure=1):
    """
    determines the equilibrium species in the given species
    input:
    input1: list of the element in the inventory
    b: inventory of the elements specified as input
    sp_g: list of the gaseous species considered
    INSERT: initial list of the condensed species, from which the iteration
    starts. This speeds up the convergence if the several equilibrium species
    are known before hand.
    total_sp_c: list of the all condensed species considered for the
    equilibrium calculation
    a: condensed part of the stoichiometric matrix
    a_g: gaseous part of the stoichiometric matrix
    trace: min. amount of the allowed mole number
    dict_of_all_sp_grt: chemical potential dictionary of the all chemical
    species at the specified temperature.
    initial_sp_c: initial list of the considered condensed chemical species
    grt_dict: chemical potential dictionary of the all chemical specis at the
    specified temperature.
    stoichiometric_dict: dictionary consisting the stoimetric data for the all
    the chemical species
    temperature: specified temperature
    v: system volume
    returns
    y: Equilibrium mole number species wise
    sp_g: list of the gaseous phase species
    sp_c: list of the condensed phase species
    """
    discarded_sp_c = []
    list_of_dependent_sp = []    # set of dependent condensed species.
    outer_iteration = 0
    while True:
        print('Outer iteration no.:', outer_iteration)

        if outer_iteration == 0:
            sp_c = INSERT.copy()
            if len(sp_c) != 0:
                for i in sp_c:
                    del total_sp_c[total_sp_c.index(i)]

            x_g = np.ones(len(sp_g)) * 10
            x_c = np.ones(len(sp_c)) * 0.1
            y_g = np.ones(len(sp_g)) * 100
            y_c = np.ones(len(sp_c)) * 0.1
            a_c = stoichiometric_coeff_matrix_generator.make_ac(
                input1, b, sp_c, stoichiometric_dict)
        else:
            a_c = stoichiometric_coeff_matrix_generator.make_ac(
                input1, b, sp_c, stoichiometric_dict)

            if np.linalg.matrix_rank(a_c.transpose()) <\
                    a_c.transpose().shape[0]:
                _, inds = sympy.Matrix(a_c.transpose()).T.rref()
                print('Dependent species is occured, checking \
                      for dependent species')

                ds2 = set(list_of_dependent_sp)
                a_c, sp_c, total_sp_c, list_of_dependent_sp = \
                    stoichiometric_coeff_matrix_generator.test_for_dependence(
                        a_c, inds, input1, b, stoichiometric_dict, sp_c,
                        dict_of_all_sp_grt, a, pis, initial_sp_c, total_sp_c)
                ds1 = set(list_of_dependent_sp)
                if ds1 == ds2:
                    dependence_count += 1
                    print('Same dependence condensed species occured',
                          dependence_count, 'times.')
                    if dependence_count > 3:
                        dependence_count = 0
                        print('Same dependence count occured more \
                              than 3 time.')
                        outer_iteration = 0
                        INSERT = list(set(sp_c.copy()) - ds2)
                        print('The following species are being inserted:',
                              INSERT)
                        sp_c = INSERT.copy()
                        try:
                            if len(sp_c) != 0:
                                for i in sp_c:
                                    del total_sp_c[total_sp_c.index(i)]
                        except:
                            pass
                        x_g = np.ones(len(sp_g)) * 1
                        y_g = np.ones(len(sp_g)) * 1
                x_c = np.ones(len(sp_c)) * 10
                y_c = np.ones(len(sp_c)) * 100
                a_c =\
                    stoichiometric_coeff_matrix_generator.make_ac(
                        input1, b, sp_c, stoichiometric_dict)
            else:
                y_c = np.ones(len(sp_c)) * 100
                x_c = np.ones(len(sp_c)) * 100
                for i in range(len(sp_c)):
                    try:
                        y_c[i] = dict_yc[sp_c[i]]
                        x_c[i] = dict_yc[sp_c[i]]
                    except:
                        y_c[i] = 1
                        x_c[i] = 1
                x_g = y_g.copy()
        if switch == 'TV':
            sp_g, sp_c, y_g, y_c, pis, a_c = sd_tv(trace, outer_iteration,
                                                   sp_c, sp_g,
                                                   x_c, x_g, y_c, y_g,
                                                   a_c, a_g, b, grt_dict, 
                                                   temperature, v)
        elif switch == 'TP':
            sp_g, sp_c, y_g, y_c, pis, a_c = sd_tp(trace, outer_iteration, 
                                                   sp_c, sp_g, x_c, x_g, y_c, 
                                                   y_g, a_c, a_g, b, grt_dict, 
                                                   temperature, v, pressure)
        else:
            print('Please enter valid switch value')

        dict_yc = {}
        for i in range(len(sp_c)):
            dict_yc[sp_c[i]] = y_c[i]

        list_of_potential = []
        for i in range(len(total_sp_c)):
            list_of_potential.append(dict_of_all_sp_grt[total_sp_c[i]] -
                                     sum(pis *
                                     a[:, initial_sp_c.index(total_sp_c[i])]))
        if (np.array(list_of_potential) > 0).all():
            break
        else:
            considered_condensed_sp = \
                (total_sp_c[list_of_potential.index(min(list_of_potential))])
            sp_c.append(considered_condensed_sp)
            del total_sp_c[list_of_potential.index(min(list_of_potential))]

        k13 = 0
        for i in range(len(y_c)):
            if y_c[i] < 0:
                discarded_sp_c.append(sp_c[k13])
                total_sp_c.append(sp_c[k13])
                del sp_c[k13]
                k13 = k13 - 1
            k13 = k13 + 1

        print('Considered condensed species at the end of one outer iteration',
              sp_c)

        outer_iteration = outer_iteration + 1

    y = []

    for i in range(len(sp_g)):
        y.append(y_g[i])
    for i in range(len(sp_c)):
        y.append(y_c[i])
    y = np.array(y)
    return y, sp_g, sp_c


def sd_tv(trace, outer_iteration, sp_c, sp_g, x_c, x_g, y_c, y_g, a_c, a_g,
          b, grt_dict, temperature, v):
    """
    input:
    trace: The species less number mole than trace will be assigned the value
    of trace, this will reduce the time require for the convergence.
    outer_iteration: Specifies the outer iteration number
    sp_c: list of the species in the condensed phase considered for the
    iteration
    sp_g: list of the species in the gaseous phase, considered for the
    iteration.
    x_c: actual mole numbers for the condensed species
    x_g: actual mole numbers for the gaseous phase species
    y_c: guessed mole numbers for the condensed species
    y_g: guessed mole numbers for the gaseous phase species
    a_c: stoichiometric matrix for the condensed species
    a_g: stoichiometric matrix for the gaseous phase species
    b  : inventory of the element as per specified in the input
    grt_dict: chemical potential dictionary
    temperature: specified temperature
    v: volume of the system
    returns:
    y_g: Equilibrium mole number gaseous species wise 
    y_c: Equilibrium mole number condensed species wise 
    sp_g: list of the gaseous phase species
    sp_c: list of the condensed phase species
    a_c: condensed part of the stoichiometric coeff matrix
    pis: values of the pi for each element
    """
    R = 8.31445984848484848484
# =============================================================================
#                           While loop
# =============================================================================
    # while iteration < 5000:
    fi = np.zeros(len(sp_g))
    old_pis = np.zeros(len(b))
    inner_iteration = 0
    while True:
        # =====================================================================
        #                        Matrix generation
        # =====================================================================
        for i in range(len(sp_g)):
            fi[i] = y_g[i] * (grt_dict[sp_g[i]] +
                              np.log(sum(y_g) * R * temperature/v/1e5) +
                              np.log(y_g[i]/sum(y_g)))    # added one

        rw1 = np.hstack((a_c.transpose(), np.zeros([len(sp_c), len(sp_c)])))

        temp_array1 = np.zeros([1, len(b)])

        for i in range(len(b)):
            sum3 = 0
            for j in range(len(y_g)):
                sum3 = sum3 + y_g[j] * a_g.transpose()[j, i]
            temp_array1[0, i] = sum3

        rw2 = np.hstack((temp_array1, np.zeros([1, len(sp_c)])))

        # r matrix preparation
        r = np.zeros([len(b), len(b)])
        for i in range(len(b)):
            for j in range(len(b)):
                r[i, j] = sum(a_g.transpose()[:, i] *
                              a_g.transpose()[:, j] * y_g)

        temp_array2 = np.zeros([len(b), len(x_c)])

        for i in range(len(b)):
            for j in range(len(x_c)):
                temp_array2[i, j] = a_c[i, j]

        rw3 = np.hstack((r, temp_array2))

        pi_matrix = np.vstack((rw1, rw3))
        temp_array4 = np.zeros([len(y_c), 1])
        for i in range(len(sp_c)):
            temp_array4[i] = grt_dict[sp_c[i]]

        temp_array6 = np.zeros([len(b), 1])
        correction = np.zeros(len(b))

        for i in range(len(b)):
            correction[i] = sum(a_g.transpose()[:, i] * y_g) - b[i]

        for i in range(len(b)):
            temp_array6[i, 0] = sum(a_g.transpose()[:, i] * fi) -\
                                correction[i]

        rhs = np.vstack((temp_array4, temp_array6))
    # =========================================================================
    #                           matrix inversion
    # =========================================================================
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                sol = solve(pi_matrix, rhs)
            except scipy.linalg.LinAlgError as err:
                if 'Singular matrix' in str(err):
                    # error handling block
                    print('Sinularity occured, breaking inner iteration loop')
                    return sp_g, sp_c, y_g, y_c, pis, a_c
            except Warning as e:     # scipy.linalg.LinAlgWarning as war
                print('Ill-conditioned matrix, breaking inner iteration')
                return sp_g, sp_c, y_g, y_c, old_pis, a_c

        pis = sol[0:len(b), 0]

        k1 = 0
        for i in range(len(b), len(b) + len(x_c)):
            x_c[k1] = sol[i, 0]
            k1 = k1 + 1

        for i in range(len(x_g)):
            x_g[i] = -fi[i] + (y_g[i] *
                               (sum(pis * a_g.transpose()[i, :]) + 1))
    # =========================================================================

    # =========================================================================
    #                           Convergence test
    # =========================================================================
        if outer_iteration == 0:
            if (abs(b - a_g.dot(y_g) - a_c.dot(y_c)) <= max(b) * 1e-5).all():
                print('break due to convergence criteria')
                break
        else:
            if abs((sum(x_g) + sum(x_c))/(sum(y_g) + sum(y_c)) - 1) < 1e-8\
                    and (abs((old_pis - pis)/pis) < 0.001).all():
                break
        old_pis = pis.copy()

    # =========================================================================
    #                           lembda determination
    # =========================================================================

        lg = np.ones(len(x_g)) * 1e15
        lc = np.ones(len(x_c)) * 1e15

        if (x_g <= 0).any():
            for i in range(len(x_g)):
                if x_g[i] <= 0:
                    lg[i] = y_g[i] / (y_g[i] - x_g[i])

        else:
            for i in range(len(x_g)):
                if x_g[i] <= 0:
                    lg[i] = y_g[i] / (y_g[i] - x_g[i])

            for i in range(len(x_c)):
                if x_c[i] <= 0:
                    lc[i] = y_c[i] / (y_c[i] - x_c[i])

        try:
            ld = np.clip(min(min(lg), min(lc)), 0, np.inf)
        except:
            ld = np.clip(min(lg), 0, np.inf)

        lembda = 0.999 * ld * (1 - ld * 0.5)
        if lembda == 0:
            lembda = 1e-3
        for i in range(len(y_g)):
            y_g[i] = y_g[i] + lembda * (x_g[i] - y_g[i])
            if abs(y_g[i]) < trace:
                y_g[i] = 1e-20
                x_g[i] = 1e-20

        for i in range(len(y_c)):

            y_c[i] = y_c[i] + lembda * (x_c[i] - y_c[i])
            if abs(y_c[i]) < trace:
                y_c[i] = 1e-20
                x_c[i] = 1e-20

        if inner_iteration > 10000:
            print('Breaking iterative loop, inner iterations exceeds \
                   iteration limit')
            break
        inner_iteration = inner_iteration + 1
    return sp_g, sp_c, y_g, y_c, pis, a_c


def sd_tp(trace, outer_iteration, sp_c, sp_g, x_c, x_g, y_c, y_g, a_c, a_g,
       b, grt_dict, temperature, v, pressure):
    """
    input:
    trace: The species less number mole than trace will be assigned the value
    of trace, this will reduce the time require for the convergence.
    outer_iteration: Specifies the outer iteration number
    sp_c: list of the species in the condensed phase considered for the
    iteration
    sp_g: list of the species in the gaseous phase, considered for the
    iteration.
    x_c: actual mole numbers for the condensed species
    x_g: actual mole numbers for the gaseous phase species
    y_c: guessed mole numbers for the condensed species
    y_g: guessed mole numbers for the gaseous phase species
    a_c: stoichiometric matrix for the condensed species
    a_g: stoichiometric matrix for the gaseous phase species
    b  : inventory of the element as per specified in the input
    grt_dict: chemical potential dictionary
    temperature: specified temperature
    pressure: pressure of the system
    returns:
    y_g: Equilibrium mole number gaseous species wise 
    y_c: Equilibrium mole number condensed species wise 
    sp_g: list of the gaseous phase species
    sp_c: list of the condensed phase species
    a_c: condensed part of the stoichiometric coeff matrix
    pis: values of the pi for each element
    """
    R = 8.31445984848484848484
    # =========================================================================
    #                           While loop
    # =========================================================================
    # while iteration < 5000:
    fi = np.zeros(len(sp_g))
    old_pis = np.zeros(len(b))
    inner_iteration = 0
    while True:
        # print('inner iteration no', inner_iteration)
        # =====================================================================
        #                        Matrix generation
        # =====================================================================
        for i in range(len(sp_g)):
            fi[i] = y_g[i] * (grt_dict[sp_g[i]] +
                              np.log(pressure) +
                              np.log(y_g[i]/sum(y_g)))    # added one

        rw1 = np.hstack((a_c.transpose(), np.zeros([len(sp_c), len(sp_c)]),
                        np.zeros([len(sp_c), 1])))

        temp_array1 = np.zeros([1, len(b)])

        for i in range(len(b)):
            sum3 = 0
            for j in range(len(y_g)):
                sum3 = sum3 + y_g[j] * a_g.transpose()[j, i]
            temp_array1[0, i] = sum3

        rw2 = np.hstack((temp_array1, np.zeros([1, len(sp_c)]),
                         np.zeros([1, 1])))

        # r matrix preparation
        r = np.zeros([len(b), len(b)])
        for i in range(len(b)):
            for j in range(len(b)):
                r[i, j] = sum(a_g.transpose()[:, i] *
                              a_g.transpose()[:, j] * y_g)

        temp_array2 = np.zeros([len(b), len(x_c)])

        for i in range(len(b)):
            for j in range(len(x_c)):
                temp_array2[i, j] = a_c[i, j]

        temp_array3 = np.zeros([len(b), 1])

        for i in range(len(b)):
            temp_array3[i, 0] = sum(a_g.transpose()[:, i] * y_g)

        rw3 = np.hstack((r, temp_array2, temp_array3))

        pi_matrix = np.vstack((rw1, rw2, rw3))
        temp_array4 = np.zeros([len(y_c), 1])
        for i in range(len(sp_c)):
            temp_array4[i] = grt_dict[sp_c[i]]

        temp_array5 = np.array([[sum(fi)]])

        temp_array6 = np.zeros([len(b), 1])
        correction = np.zeros(len(b))

        for i in range(len(b)):
            correction[i] = sum(a_g.transpose()[:, i] * y_g) - b[i]

        for i in range(len(b)):
            temp_array6[i, 0] = sum(a_g.transpose()[:, i] * fi) -\
                                correction[i]

        rhs = np.vstack((temp_array4, temp_array5, temp_array6))
        # print(sp_g, rhs, pi_matrix)
    # =========================================================================
    #                           matrix inversion
    # =========================================================================
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                # sol = np.linalg.inv(pi_matrix).dot(rhs)
                sol = solve(pi_matrix, rhs)
            except scipy.linalg.LinAlgError as err:
                if 'Singular matrix' in str(err):
                    # your error handling block
                    print('sinularity occured, breaking inner iteration')
                    return sp_g, sp_c, y_g, y_c, pis, a_c
            except Warning as e:     #  scipy.linalg.LinAlgWarning as war
                print('Ill-conditioned matrix, breaking inner iteration', e)
                return sp_g, sp_c, y_g, y_c, old_pis, a_c

        pis = sol[0:len(b), 0]

        k1 = 0
        for i in range(len(b), len(b) + len(x_c)):
            x_c[k1] = sol[i, 0]
            k1 = k1 + 1

        u = sol[-1]

        # print('u:', u)

        for i in range(len(x_g)):
            x_g[i] = -fi[i] + (y_g[i] *
                               (u + 1 + sum(pis * a_g.transpose()[i, :])))
        # y_c = x_c.copy()
    # =========================================================================

    # =========================================================================
    #                           Convergence test
    # =========================================================================
        if outer_iteration == 0:
            if (abs(b - a_g.dot(y_g)) - a_c.dot(y_c) <= max(b) * 1e-6).all():
                print('break due to convergence criteria')
                break
        else:
            # if (abs((old_pis - pis)/pis) < 0.001).all():
            #     print('pie convergence criteria met')
            #     break
            if abs(u) < 1e-8 and \
            abs((sum(x_g) + sum(x_c))/(sum(y_g) + sum(y_c)) - 1) < 1e-8:
                # break
                break
        old_pis = pis.copy()


        # rmse = np.sqrt(np.sum((y_g - x_g))**2 / len(x_g))
        # print(rmse)
        # if rmse < 1e-2:
        #     break
    # =========================================================================
    #                           lembda determination
    # =========================================================================

        lg = np.ones(len(x_g)) * 1e15
        lc = np.ones(len(x_c)) * 1e15

        if (x_g <= 0).any():
            for i in range(len(x_g)):
                if x_g[i] <= 0:
                    lg[i] = y_g[i] / (y_g[i] - x_g[i])
                    # if lg[i] <= 0:
                    #     lg[i] = 1e10
        else:
            for i in range(len(x_g)):
                if x_g[i] <= 0:
                    lg[i] = y_g[i] / (y_g[i] - x_g[i])
                    # if lg[i] <= 0:
                    #     lg[i] = 1e10
            for i in range(len(x_c)):
                if x_c[i] <= 0:
                    lc[i] = y_c[i] / (y_c[i] - x_c[i])
                    # if lc[i] <= 0:
                    #     lc[i] = 1e10

        try:
            ld = np.clip(min(min(lg), min(lc)), 0, np.inf)
        except:
            ld = np.clip(min(lg), 0, np.inf)

        lembda = 0.999 * ld * (1 - ld * 0.5)

        if lembda == 0 or abs(lembda) >= 1e10:
            lembda = 1e-3
            # print('hi', lembda)
        # print(lembda)
        # if inner_iter == 0:
        for i in range(len(y_g)):
            y_g[i] = y_g[i] + lembda * (x_g[i] - y_g[i])
            if abs(y_g[i]) < trace:
                y_g[i] = 1e-20
                x_g[i] = 1e-20

        for i in range(len(y_c)):

            y_c[i] = y_c[i] + lembda * (x_c[i] - y_c[i])
            if abs(y_c[i]) < trace:
                y_c[i] = 1e-20
                x_c[i] = 1e-20

        # print(list(zip(sp_g, y_g)), list(zip(sp_c, y_c)))
        if inner_iteration > 15000:
            print('program break from iteration limit')
            break
        # print(y_g)
        inner_iteration = inner_iteration + 1
    return sp_g, sp_c, y_g, y_c, pis, a_c


def min_fun_helmholtz(x, species, grt_dict, temperature, v):
    """
    helmholtz function
    input:
    x: array containing mole numbers
    species: list of species
    grt_dict: dictionary of chemical potential for all the chemical species
    temperature: specified temperature
    v: system volume
    returns: helmholtz function value
    """
    R = 8.31445984848484848484   # gas constant
    sum1 = 0
    gas = 0
    liquid = 0
    cr = 0
    condense_phase_match = '\([A|B|H|L|a|b|c|d|d\'|e|X|cr|I|III|II]*\)'
    for i in range(len(species)):

        # finds total moles of gaseous species, liquid and crystal form
        if re.findall(condense_phase_match, species[i]) == []:
            gas = gas + x[i]

        elif re.findall(condense_phase_match,
                        species[i]) == ['(L)']:
            liquid = liquid + x[i]

        else:
            cr = cr + x[i]

    for i in range(len(species)):
        if re.findall(condense_phase_match, species[i]) == []:
            sum1 = sum1 + x[i] * (grt_dict[species[i]]
                                  + np.log(gas * R * 1e-5 *
                                           temperature / v)
                                  + np.log(x[i] / gas))
        else:
            sum1 = sum1 + x[i] * (grt_dict[species[i]])
    return sum1


def gibbs_calculate(x, species, grt_dict, temperature, P):
    """
    Gibbs function
    input:
    x: array containing mole numbers
    species: list of species
    grt_dict: dictionary of chemical potential for all the chemical species
    temperature: specified temperature
    P: system pressure
    returns: Gibbs function value
    """
    sum1 = 0
    gas = 0
    liquid = 0
    cr = 0
    condense_phase_match = '\([A|B|H|L|a|b|c|d|d\'|e|X|cr|I|III|II]*\)'
    for i in range(len(species)):

        # finds total moles of gaseous species, liquid and crystal form
        if re.findall(condense_phase_match, species[i]) == []:
            gas = gas + x[i]

        elif re.findall(condense_phase_match,
                        species[i]) == ['(L)']:
            liquid = liquid + x[i]

        else:
            cr = cr + x[i]

    for i in range(len(species)):
        if re.findall(condense_phase_match, species[i]) == []:
            sum1 = sum1 + x[i] * (grt_dict[species[i]]
                                  + np.log(P)
                                  + np.log(x[i] / gas))
        else:
            sum1 = sum1 + x[i] * (grt_dict[species[i]])

    return sum1


def jac_eq_const(x):
    return a
