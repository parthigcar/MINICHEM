================
Running the code
================

Before running the code, please check :ref:`Requirements`, for the necessary python packages required to run the code. To run this program following files are required:

- chem_parse.py
- data_process.py
- mini.py
- plot_hv.py
- pythormoread.py
- rf.py
- species_search.py
- stoichiometric_coeff_matrix_generator.py

MINICHEM uses following input file.

**Input file:**

- thermo_chemical_database.txt


The code generates following output files.

**Output files:**

- iom.txt: Gives the released mole, RF, and details about considered and released species.
- released_sp.txt: Gives the released species mole
- released_mole_el.txt: element wise result of the released mole.


------------------------
Defining input inventory
------------------------

The input inventory can be specified by modifying the input1 list in ``minichemp.py``. Each list input element should be specified as string starting with its the mole number followed by the element name.

For example,

.. code-block::
    :linenos:

    input1 = ['2H', '1N', '1O']

-------------------------------------------------------
Defining temperature, pressure, volume parameters
-------------------------------------------------------

For each case ((T, P) or (T, V)), the temperature (in K), P (in bars) and V (in :math:`m^3`) needs to specified, even if one of these parameter is redundant.

------
Method
------

MINICHEM can calculate (T, P) or (T, V) chemical equilibrium using two method:
1. Quadratic gradient descent minimising method
2. Sequential Quadratic Programming

For this, there are two switch provided named ``method`` and ``switch``
- The valid value for the ``method`` is ``SLSQP`` or ``MINI`` (which calculates chemical equilibrium using `Quadratic gradient descent minimising method`)
- The valid value for the ``switch`` is ``TP`` or ``TV``.

-----
Trace
-----

By default the value of ``trace`` is set to :math:`10^{-25}`.