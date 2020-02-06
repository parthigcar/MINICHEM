------------------
Extending database
------------------

Since, the thermochemical database is fed via dictionary (``dict_of_all_sp_grt``), the thermochemical database can be easily modified by either providing NASA 9 polynomial in ``thermo_dict`` or directly providing chemical potentials to ``dict_of_all_sp_grt`` at specified temperature.


.. code-block::
    :linenos:

    thermo_dict, stoichiometric_dict = pythermoread.thermoread()
    dict_of_all_sp_grt = pythermoread.calculate_grt(grt_dict, temperature, thermo_dict)

----------
References
----------
.. [Kenneth1956] Kenneth Denbigh, 1956. The Principles of Chemical Equilibrium, with Applications in Chemistry and Chemical Engineering., :math:`4^{th}` Edition. Cambridge University Press, Cambridge.
.. [Bertsekas2014] Bertsekas, D. P., 2014. Constrained Optimization and Lagrange Multiplier Methods. Academic Press.
.. [Weisstein2020] Weisstein, E. W., 2020. Lagrange Multiplier. http://mathworld.wolfram.com/LagrangeMultiplier.html.
.. [White1958] White, W. B., Johnson, S. M., Dantzig, G. B., 1958. Chemical Equilibrium in Complex Mixtures. The Journal of Chemical Physics 28 (5), 751–755.
.. [McBride_and_Gordon1996a] McBride, B. J., Gordon, S., 1996b. Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications-II: User Manual and Program Description. Vol. 2. National Aeronautics and Space Administration, Office of Management, Scientific and Technical Information Program, Cleveland, Ohio, 44135.
.. [Eriksson1979] Eriksson, G., 1979. An algorithm for the computation of aqueous multi-component, multiphase equilibria. Analytica Chimica Acta 112 (4), 375–383.
.. [Eriksson_and_Hack1990] Eriksson, G., Hack, K., 1990. ChemSage—a computer program for the calculation of complex chemical equilibria. Metallurgical transactions B 21 (6), 1013–1023.
.. [Gunnar_Eriksson1971] Gunnar Eriksson, 1971. Thermodynamic Studies of High Temperature Equilibria. III. SOLGAS, a Computer Program for Calculating the Composition and Heat Condition of an Equilibrium Mixture. Acta Chemica Scandinavica 25, 2651–2658.