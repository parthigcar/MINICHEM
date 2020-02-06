.. MINICHEM documentation master file, created by
   sphinx-quickstart on Thu Jan  9 08:54:19 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MINICHEM (MINImisation of CHEMical potentials)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The MINICHEM (MINImisation of CHEMical potentials) is a code to calculate the In-vessel source term using thermochemical equilibrium approach. The code calculates the equilibrium species in the reactor vessel during ULOFA. The released inventory in the cover gas is expressed in terms of release fraction. Although, the
present code only depicts the chemical aspects of the release behaviour of RN,
which we consider it as a first step towards mechanistic model development for
oxide fuelled SFRs. For the purpose of this analysis, ULOF event resulting in
whole core melt is considered.
The in-vessel source term is determined using chemical equilibrium approach
with no mixture assumption. No mixture assumption essentially means that
during chemical equilibrium, the mixing properties of the species are not
considered; The estimated equilibrium species mole numbers corresponds to
the vapour pressure of the species at specified temperature. This assumption
usually leads to conservative estimates. With the help of this equilibrium
species, distribution of RN in three phases (solid, gas and liquid) can be
determined; From this information, release fractions of RNs to the cover gas
are evaluated. The code is capable to calculate the equilibrium species at both (T, P) and (T, V) case.

Apart from source term calculations, the code can be used for the
general purpose chemical equilibrium calculations. The
code is developed at the Reactor Shielding and Data Division (RSDD) at `Indira Gandhi Center for Atomic Research, India <http://www.igcar.gov.in/>`_. The source code is hosted on `GitHub <https://github.com/geekonloose/MINICHEM>`_.

.. admonition:: Recommended publication for citing
   :class: tip

   Patel, P.R., Arul, A.J., In-vessel source term calculation for a hypothetical core disruptive accident using chemical equilibrium approach for a medium sized sodium cooled fast reactor. Submitted to Nuclear Engg. & Design 25.


.. only:: html

   --------
   Contents
   --------

.. toctree::
    :maxdepth: 1

    Requirements
    theory
    code_structure
    running
    extending_database

.. meta::
    :google-site-verification: AVg7xNWpclGRAt67p2S1ffN28JPdG0qDyMjl1BvhFLU
