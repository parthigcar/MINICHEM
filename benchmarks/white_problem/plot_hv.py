'''
File: plot_hv.py
Project: MINICHEM
File Created: Tuesday, 7th January 2020 8:49:51 am
Author: John Arul & Parth (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Last Modified: Thursday, 9th January 2020 8:24:52 am
Modified By: John Arul & Parth Patel (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Copyright: IGCAR - 2020
'''

import numpy as np
import pythermoread
import holoviews as hv
import matplotlib.pyplot as plt
import networkx as nx
from holoviews import opts, renderer, render
from holoviews import dim
from bokeh.plotting import show


import pandas as pd
from itertools import chain
# import from other modules
# from functions import *
import re
# thermo_dict, stoichiometric_dict = pythermoread.thermoread()
# grt_dict = {}
# temperature = 900
# dict_of_all_sp_grt = pythermoread.calculate_grt(grt_dict, temperature,
# thermo_dict)
# adding lanthanides and actinides oxide and elemental thermochemical database
# h_rt_dict, s_r_dict, temp_range = make_sr_hrt_dict()
# st_dict = make_stoichiometric_coeff()
# stoichiometric_dict, dict_of_all_sp_grt = data_base(temperature,
# dict_of_all_sp_grt, stoichiometric_dict)
# for i in st_dict:
#     stoichiometric_dict[i] = st_dict[i]
# input1 = ['He', 'Kr', 'Xe', 'O', 'H', 'I', 'Br', 'Na', 'Cs', 'Rb', 'Ba',
#           'Sr', 'Mo', 'Zr', 'U', 'Sb']
# include_el = ['I', 'Cs']   # ['U', 'Mo', 'Xe', 'Zr']
# include_sp = []


def plot_hv(input1, stoichiometric_dict, include_el, opfilename, Min=0, 
            Max=1e9):
    """
    Plotting module
    input1: list of input elements
    include_el: list of element for which sankey chart is drawn
    Min: min mole number species to be included in chart
    Max: max mole number species to be included in the chart
    opfilename: opfilename with extension (allowed formats:
    ['html', 'json', 'auto', 'png', 'widgets', 'scrubber', 'auto', None]) 
    returns:
    saves sankey chart
    """
    f = open('released_sp.txt', 'r')
    data = f.readlines()
    sp_dict = {}
    for line in data:
        cols = line.split()
        sp_dict[cols[0]] = float(cols[1])
    e = []
    sp = []
    val = []
    hv.extension('bokeh')
    for i in include_el:
        for j in sp_dict:
            s1 = set()
            for j1 in stoichiometric_dict[j]:
                s1.add(j1[0])

            if i.upper() in s1:
                if i == j:
                    if sp_dict[j] > Min and sp_dict[j] < Max:
                        e.append(i)
                        sp.append(j + '*')
                        val.append(sp_dict[j])
                else:
                    if sp_dict[j] > Min and sp_dict[j] < Max:
                        e.append(i)
                        sp.append(j)
                        val.append(sp_dict[j])

    data = {'elements': e, 'species': sp, 'values': val}
    print(data)
    df = pd.DataFrame(data)
    print(df)
    sankey = hv.Sankey(df)

# , show_values=False
    sankey.opts(label_position='right', edge_color='species',
                node_color='elements',
                labelled=['elements', 'species'], xlabel='part',
                show_values=False)
    hv.save(sankey, opfilename, backend='bokeh')
    # show(hv.render(sankey))
    # hv.extension('matplotlib')
    # hv.save(sankey, 'curve1.png', fmt='svg')
    # plt.show(hv.render(sankey))


# plot_hv(input1, sp_dict, include_el)
