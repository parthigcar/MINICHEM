'''
File: pythermoread
Project: MINICHEMP
File Created: Monday, 8th April 2019 5:33:44 pm
Author: John Arul & Parth (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Last Modified: Monday, 8th April 2019 5:34:50 pm
Modified By: John Arul & Parth Patel (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Copyright: IGCAR - 2019
'''
import numpy as np
import itertools
import math
import re
from scipy.interpolate import InterpolatedUnivariateSpline
from chem_parse import chemparse
# Functions


def HRT(a1, a2, a3, a4, a5, a6, a7, b1, b2, t):
    """
    Finds the value of H/RT
    INPUT:
    Cp coefficient: a1, a2, a3, a4, a5, a6, a7,
    b1, b2: Integration coefficients
    """
    return -a1 * t**(-2) + a2 * t**(-1) * math.log(t) + a3 + a4 * t/2 \
        + a5 * t**(2) / 3 + a6 * t**(3) / 4 + a7 * t**(4) / 5 + b1 / t


def SR(a1, a2, a3, a4, a5, a6, a7, b1, b2, t):
    """
    Finds the value of S/R
    INPUT:
    Cp coefficient: a1, a2, a3, a4, a5, a6, a7,
    b1, b2        : Integration coefficients
    t             : desired temperature
    """
    return -a1*t**(-2)/2 - a2*t**(-1) + a3*math.log(t) + a4 * t + a5 * t**(2)/\
    2 + a6 * t**(3) / 3 + a7 * t**(4) / 4 + b2


def GRT1(a1, a2, a3, a4, a5, a6, a7, b1, b2, temp):
    """
    Calculates the thermochemical potential from the spcified 9 polynomial
    coefficients and the temperature information.
    Returns:
    Chemical potential at the specified temperature.
    """
    # return HRT(a1, a2, a3, a4, a5, a6, a7, b1, b2, temp) -\
    #    SR(a1, a2, a3, a4, a5, a6, a7, b1, b2, temp)
    return HRT(a1, a2, a3, a4, a5, a6, a7, b1, b2, temp) -\
        SR(a1, a2, a3, a4, a5, a6, a7, b1, b2, temp)


def thermoread():
    """
    From the make lib input, this function reads all NASA 9 polynomial
    thermochemical potentials for the all the chemical species and converts
    this database into the dictionary.
    This function also returns the stoichiometric data for all the
    thermochemical species specified in the dictionary.
    """
    f = open('thermo_chemical_database.txt', 'r')   # processed from NASA lib
    data = f.readlines()

    stoichiometric_dict = {}
    thermo_dict = {}
    for line in data:
        cols = line.split()

        # if re.search('\(cr\)|\(L\)|\(a\)|\(b\)|\(c\)', cols[0]):
        if float(cols[13]) > 0:    # phase == condensed
            sp_name = cols[0].replace('*', '')
            lis = cols[3:13]    # stoichiometric information

            for i in range(5):
                if sp_name not in stoichiometric_dict.keys():
                    stoichiometric_dict[sp_name] = []
                stoichiometric_dict[sp_name].append(
                    [lis[2 * i], float(lis[2 * i + 1])])
            tl = [float(cols[14]), float(cols[15])]
            mwt = float(cols[16])
            thermo = list(map(float, cols[17:26]))
            # thermo_dict = [[[tl1, tl2], [tl3,tl4]],
            #                [[coeff_data1], [coeff_data2]]]
            if sp_name not in thermo_dict.keys():
                thermo_dict[sp_name] = [[], []]

            thermo_dict[sp_name][0].append(tl)
            thermo_dict[sp_name][1].append(thermo)
        else:    # else phase is gas
            sp_name = cols[0].replace('*', '')
            lis = cols[3:13]    # stoichiometric information
            # print(lis)
            for i in range(5):
                if sp_name not in stoichiometric_dict.keys():
                    stoichiometric_dict[sp_name] = []
                stoichiometric_dict[sp_name].append(
                    [lis[2 * i], float(lis[2 * i + 1])])
            ntl = float(cols[1])
            # print(ntl, 'ntl')
            temp_list = []
            for i1 in range(int(ntl)):
                temp_list.append([])
                temp_list[i1].append(float(cols[14 + 2 * i1]))
                temp_list[i1].append(float(cols[14 + (2 * i1 + 1)]))
            if sp_name not in thermo_dict.keys():
                thermo_dict[sp_name] = [[], []]

            index = 14 + (2 * int(ntl) - 1) + 2
            thermo_dict[sp_name] = [temp_list,
                                    [list(map(float, cols[index:index + 9])),
                                     list(map(float,
                                              cols[index + 9: index + 18])),
                                     list(map(float,
                                              cols[index + 18: index + 27]))]]
    return thermo_dict, stoichiometric_dict


def calculate_grt(grt_dict, input_temp, thermo_dict):
    """
    The function calculates the chemical potential for the all the NASA CEA
    lib chemical specis at specified input temperature.
    input:
    grt_dict: This dictionary is empty
    input_temp: input temperature at which the chemical potential will be
    calculated.
    thermo_dict: dictionary containing NASA 9 polynomial thermochemical
                 database.
    Returns:
    Dictionary containing the chemical potentials at the specified temperature.
    """
    for sp_name in thermo_dict.keys():
        temp_list = thermo_dict[sp_name][0]
        # if len(temp_list)>3:
        #     print(sp_name, len(temp_list))
        temp_range = \
            list(itertools.chain.from_iterable(thermo_dict[sp_name][0]))
        k1 = 0
        if min(temp_range) <= input_temp <= max(temp_range):
            for i1 in temp_list:

                if i1[0] <= input_temp <= i1[1]:
                    index1 = k1
                    break
                k1 = k1 + 1
            a1, a2, a3, a4, a5, a6, a7, b1, b2 =\
                thermo_dict[sp_name][1][index1]
            grt_dict[sp_name] = GRT1(a1, a2, a3, a4, a5,
                                     a6, a7, b1, b2, input_temp)
        else:
            # This portion of code needs to be improvised
            # (extrapolation of the condensed species)
            # print('temp_range is outside of given ranges',
            #       sp_name, min(temp_range), max(temp_range))
            if input_temp < min(temp_range):
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = thermo_dict[sp_name][1][0]
                grt_dict[sp_name] = GRT1(a1, a2, a3, a4,
                                         a5, a6, a7, b1, b2, min(temp_range))
                # grt_dict[sp_name] = GRT1(a1, a2, a3, a4, a5,
                #                          a6, a7, b1, b2, input_temp)
            else:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 =\
                    thermo_dict[sp_name][1][-1]
                # grt_dict[sp_name] = GRT1(a1, a2, a3, a4, a5,
                #                          a6, a7, b1, b2, max(temp_range))
                # grt_dict[sp_name] = GRT1(a1, a2, a3, a4, a5,
                #                          a6, a7, b1, b2, input_temp)
                # f =\
                # InterpolatedUnivariateSpline([max(temp_range), input_temp],
                #                              [GRT1(a1, a2, a3, a4, a5, a6,
                #                                    a7, b1, b2,
                #                                    max(temp_range)),
                #                               GRT1(a1, a2, a3, a4, a5, a6,
                #                                    a7, b1, b2, input_temp)],
                #                               k=1)
                # This is hack for condensed species
                grt_dict[sp_name] = 1e6       # f(input_temp)
    return grt_dict


def grt(g_hts, hf, T):
    """
    Function to calculate the chemical potentials from the janaf table values,
    and prints the chemical potential for that temperature.
    """
    return print(-g_hts/R + hf/R/T * 1000)


def only_grt(grt_dict, strlist):
    '''
    takes the input of the complete combination of the input element as
    dictionary and the list of desired elements which we want to keep in
    dictionary. The function will delete other species combination.
    Input:
    grt_dict: all combination of input1 from thermochem lib
    strlist: list of the desired elements
    Output:
    grt_dict: only g/rt data of the desired elements which are in strlist.
    '''
    for i in list(grt_dict):
        if i not in strlist:
            del grt_dict[i]
    return grt_dict
