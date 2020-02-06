.. MINICHEM documentation master file, created by
   sphinx-quickstart on Thu Jan  9 08:54:19 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MINICHEM (MINImisation of CHEMical potentials)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The MINICHEM (MINImisation of CHEMical potentials) is a code to calculate the  thermochemical equilibrium of various elements present in a nuclear reactor core. The code calculates the equilibrium species in a given volume and temperature or pressure and temperature. It can treat solid, liquid and gas phases. Phases in a solid are not treated. The program has been used to calculate inventory released to the cover gas from a molten core of a sodium cooled fast reactor. 

The substances assumed to be in a no-mixture or real mixture state. 
No mixture assumption essentially means that
during chemical equilibrium, the mixing properties of the species are not
considered; The estimated equilibrium species mole numbers corresponds to
the vapour pressure of the species at specified temperature. This assumption
usually leads to conservative estimates. With the help of this equilibrium
species, distribution of RN in three phases (solid, gas and liquid) can be
determined; From this information, release fractions of RNs to the cover gas
are evaluated. To invoke the real mixture assumption, the excess functions of the corresponding species are required in the database.

Apart from source term calculations, the code can be used for the
general purpose chemical equilibrium calculations with appropriate database. The project was started under the Reactor Shielding and Data Division (RSDD) at `Indira Gandhi Center for Atomic Research, India <http://www.igcar.gov.in/>`_ and an IAEA CRP on the subject. The source code is hosted on `GitHub <https://github.com/geekonloose/MINICHEM>`_. MINICHEM is distributed under GNU General Public License v3.0 

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
