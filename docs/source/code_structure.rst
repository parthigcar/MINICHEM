==============
Code structure
==============
The MINICHEM is written in Python 3.7.1. However it is also tested with
Python >= 3.7.1. Main modules of the MINICHEM are given in the
:numref:`code-structure`. Each of module is explained in the subsequent section.


.. _code-structure:

.. figure:: ../_images/code_structure_low_res.png
   :align: center


-------
Modules
-------

^^^^^^^^^
PyMakelib
^^^^^^^^^

``Filename: pymakelib.py``

This module reads ``NASA's CEA thermochemical database`` file named ``thermo.inp``. This module is translated version of the original Fortran ``makelib.f`` module (available in CEA code package). The data for the gas phase is extrapolated for the higher temperature. The restructured data is stored in ``thermochemical_database.txt``.

^^^^^^^^^^^^
PyThermoread
^^^^^^^^^^^^

``Filename: pythermoread.py``

This module reads extrapolated thermochemical databases in ``thermochemical_database.txt`` and calculates the chemical potentials and stoichiometric values at the specified temperature. These calculated chemical potentials for the respective species are stored in the dictionary named as: ``grt_dict`` and ``stoichiometric_dict``. Following are the functions are part of the pythermoread.

.. function:: HRT(a1, a2, a3, a4, a5, a6, a7, b1, b2, t):

    Finds the value of :math:`\frac{H}{RT}`

    :param a1:  :math:`C_p` coefficient
    :param a2:  :math:`C_p` coefficient
    :param a3:  :math:`C_p` coefficient
    :param a4:  :math:`C_p` coefficient
    :param a5:  :math:`C_p` coefficient
    :param a6:  :math:`C_p` coefficient
    :param a7:  :math:`C_p` coefficient
    :param b1: Integration coefficient
    :param b2: Integration coefficient
    :param t: Temperature (K)
    :returns: Returns the value of :math:`\frac{H}{RT}`


.. function:: SR(a1, a2, a3, a4, a5, a6, a7, b1, b2, t):

    Finds the value of :math:`\frac{S}{R}`

    :param a1:  :math:`C_p` coefficient
    :param a2:  :math:`C_p` coefficient
    :param a3:  :math:`C_p` coefficient
    :param a4:  :math:`C_p` coefficient
    :param a5:  :math:`C_p` coefficient
    :param a6:  :math:`C_p` coefficient
    :param a7:  :math:`C_p` coefficient
    :param b1:  Integration coefficient
    :param b2:  Integration coefficient
    :param t:   Temperature (K)
    :returns: Returns the value of :math:`\frac{S}{R}`


.. function:: GRT1(a1, a2, a3, a4, a5, a6, a7, b1, b2, temp):

    Calculates the thermochemical potential from the specified 9 polynomial
    coefficients and the temperature information.

    :param a1:  :math:`C_p` coefficient
    :param a2:  :math:`C_p` coefficient
    :param a3:  :math:`C_p` coefficient
    :param a4:  :math:`C_p` coefficient
    :param a5:  :math:`C_p` coefficient
    :param a6:  :math:`C_p` coefficient
    :param a7:  :math:`C_p` coefficient
    :param b1:  Integration coefficient
    :param b2:  Integration coefficient
    :param temp:   Temperature (K)
    :returns: Returns thermochemical potential from the specified 9 polynomial coefficient and temperature.


.. function:: thermoread():

    From the ``thermochemical_database.txt``, this function reads all NASA 9 polynomial thermochemical potentials for the all the chemical species and converts this database into the dictionary. This function also returns the stoichiometric data for all the thermochemical species.

    :returns: ``thermo_dict``, ``stiochemitric_dict``. Dictionary containing all NASA 9 polynomial coefficients and stoichiometric coefficient information for all chemical species specified in ``thermochemical_database.txt``.


.. function:: calculate_grt(grt_dict, input_temp, thermo_dict):

    The function calculates the chemical potential using ``thermo_dict`` at specified input temperature and returns in the form of dictionary.

    :param grt_dict: Dictionary to store the thermochemical potentials at specified temperature
    :param input_temp: input temperature at which the chemical potential to be calculated.
    :param thermo_dict: dictionary containing NASA 9 polynomial thermochemical
                 database.
    :returns: ``grt_dict``. Dictionary containing the chemical potentials at the specified temperature.


.. function:: only_grt(grt_dict, strlist):

    This function provides functionality to calculate the chemical equilibrium for the desired chemical species only. The function takes the input of the complete combination of the input element as dictionary and the list of desired species which we want to calculate the thermochemical equilibrium. The function will delete other species combination.

    :param grt_dict: all combination of input1 from thermochem lib
    :param strlist: list of the desired elements
    :returns: ``grt_dict``, updated ``grt_dict``, which only contains :math:`\frac{g}{RT}` data of the desired elements which are in ``strlist``.

^^^^^^^^^^^^^^
Species search
^^^^^^^^^^^^^^

``Filename: species_search.py``

In order to calculate the chemical equilibrium using the given input chemical elements/species, list of possible species which are combination of the input chemical elements/species. This module finds the combination of the input chemical elements/species from the grt_dict species list.


.. function:: combination_search(species, grt_dict, combination_sp):

    Searches the combination of the input elements in the ``grt_dict``

    :param species: List of species for which combination search will be taken out
    :param grt_dict: Dictionary of the chemical potentials
    :param combination_sp: list of combination of species
    :returns: ``combination_sp``, list of species containing combination of species list.
    :Example:

        If species is H2O,
        Then the combinations in ``grt_dict`` might be: OH, H2O2, H, O2 etc.

.. function:: el(species):

    The function takes the input species (list form), and convert
    the element of the list which can be compound/element to the element.
    INPUT: species (list)

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Stoichiometric coefficient matrix generator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Filename: stoichiometric_coeff_matrix_generator.py``

This module generates the stoichiometric coefficient. The module contains following functions:

.. function:: stoi(species, input1, stoichiometric_dict):

    Makes one row of the stoichiometry coefficient.

    :param species:List of species (species containing combination of the input elements)
    :param input1: list of the elements provided as input.
    :param stoichiometric_dict: dictionary of the species with the information
    :returns:list of the row of the stoichiometric matrix


.. function:: make_ac(input1, b, considered_sp_c, stoichiometric_dict):

    Makes stoichiometry matrix for the condensed species.

    :param input1: list of the input elements (provided by user)
    :param b: input element inventory value (provided by user)
    :param considered_sp_c: set of considered species in the condensed phase
    :returns a_c: condensed stoichiometry matrix

.. function:: test_for_dependence(a_c, inds, input1, b, stoichiometric_dict,
                        sp_c, dict_of_all_sp_grt, a, pis, initial_sp_c,
                        total_sp_c):

    Sometimes the two or more dependent species in the condensed stoichiometric coefficient matrix can occurs, this can make the condensed stoichiometric coefficient matrix singular. This function checks the ``a_c`` matrix for the dependent row, and returns the list of the
    dependent rows as well as the list of the dependent species.

    :param a_c: condensed part of stoichiometric coefficient matrix
    :param inds: Indices which are found to be dependent (calculated from the reduced row echelon form)
    :param input1: list of the input elements (provided by user)
    :param b: input element inventory value (provided by user)
    :param stoichiometric_dict: dictionary of the species with the information
    :param sp_c: list of the condensed chemical species
    :param dict_of_all_sp_grt: dictionary containing chemical potential of the all chemical species at the specified temperature
    :param a: gas part of the stoichiometric coefficient matrix
    :param pis: list of the :math:`\pi _i`
    :param initial_sp_c: list of the all condensed species in the dict_of_all_sp_grt
    :param total_sp_c: list of the all condensed chemical species which are being considered for the equilibrium calculation

    :returns: updated ``a_c``, updated list ``sp_c``, updated list ``total_sp_c``, returns the list of dependent species in the a_c matrix

^^^^
MINI
^^^^

``Filename: mini.py``

This module contains necessary functions to calculate the thermochemical equilibrium using the ``Quadratic gradient descent minimisation method`` and ``SLSQP`` method. The calculation using SLSQP method is performed using built in ``scipy`` module named ``scipy.optimize.optimize``. The major functions ``MINI`` module are described below:

.. function:: mini_solver(input1, b, sp_g, INSERT, total_sp_c, a, a_g, trace,
                          dict_of_all_sp_grt, initial_sp_c, grt_dict,
                          stoichiometric_dict, switch,
                          temperature, v=0, pressure=1):

    Determines the equilibrium species in from given input elment/species list. The function contain ``sd_tv`` and ``sd_tp`` sub-functions, which basically calculates the equilibrium species for the (T, V) and (T, P) cases respectively.

    :param input1: list of the element in the inventory
    :param b: inventory of the elements specified as input
    :param sp_g: list of the gaseous species considered
    :param INSERT: initial list of the condensed species, from which the iteration starts. This speeds up the convergence if the several equilibrium species are known before hand.
    :param total_sp_c: list of the all condensed species
     considered for the equilibrium calculation
    :param a: condensed part of the stoichiometric matrix
    :param a_g: gaseous part of the stoichiometric matrix
    :param trace: min. amount of the allowed mole number
    :param dict_of_all_sp_grt: chemical potential dictionary of the all chemical species at the specified temperature.
    :param initial_sp_c: initial list of the considered condensed chemical species
    :param grt_dict: chemical potential dictionary of the all chemical specis at the specified temperature.
    :param stoichiometric_dict: dictionary consisting the stoimetric data for the all the chemical species
    :param temperature: specified temperature
    :param v: system volume
    :returns: ``y``, Equilibrium mole number species wise. ``sp_g``, list of the gaseous phase species. ``sp_c``, list of the condensed phase species.


.. function:: min_fun_helmholtz(x, species, grt_dict, temperature, v):

    This function calculates helmholtz function for the guessed array x containing mole numbers at each iteration

    :param x: array containing mole numbers
    :param species: list of species
    :param grt_dict: dictionary of chemical potential for all the chemical species
    :param temperature: specified temperature
    :param v: system volume
    :returns: helmholtz function value for x


.. function:: gibbs_calculate(x, species, grt_dict, temperature, P):

    This function calculates Gibbs function for the guessed array x containing mole numbers at each iteration

    :param x: array containing mole numbers
    :param species: list of species
    :param grt_dict: dictionary of chemical potential for all the chemical species
    :param temperature: specified temperature
    :param P: system pressure
    :returns: Gibbs function value

^^
RF
^^

``Filename: rf.py``

This module takes the output equilibrium mole number array as input and calculates the release fractions in the cover gas.


.. function:: rf(y, species, input1, stoichiometric_dict, el_inventory):

    Takes the output mole number array and returns the dictionary with the
    release fraction in the cover gas.

    :param y: output array containing mole number
    :param species: list of the species considered
    :param input1: list of the element initially considered.
    :param stoichiometric_dict: dictionary containing the stoichiometric information for the all the chemical species.
    :param el_inventory: dictionary containing the information
     about the input inventory specified.
    :returns: Prints the cover gas release fractions and writes the output in
     the iom.txt, released_mole_el.txt, released_sp.txt

^^^^^^^^
Plotting
^^^^^^^^

``Filename: plot_hv.py``

This module plots the sankey charts using ``holoviews`` module.

.. function:: plot_hv(input1, stoichiometric_dict, include_el, opfilename,
                      Min=0, Max=1e9):

    Plotting module

    :param input1: list of input elements
    :param include_el: list of element for which sankey chart is drawn
    :param Min: min mole number species to be included in chart
    :param Max: max mole number species to be included in chart
    :param opfilename: opfilename
    :returns: saves sankey chart