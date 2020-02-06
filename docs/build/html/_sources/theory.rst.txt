------
Theory
------

The problem of determining the equilibrium concentration of various species
given a set of input elements can be formulated in terms of either entropy,
Gibbs or Helmholtz function. For example, if the system is defined in terms of
temperature and pressure, then minimisation of Gibbs function is appropriate
since temperature and pressure are independent variables of Gibbs function
[Kenneth1956]_ . If the system is defined in terms of
temperature and volume, then minimisation of Helmholtz function is appropriate.

For a multi-phase system with :math:`\mathrm{(x_1, x_2,\ldots, x_{N_g},\ldots
x_{N_g + N_s})}` moles of :math:`\mathrm{(N_g + N_s)}` species containing
:math:`\mathrm{N_e}` elements, the Helmholtz function to be minimized is as
follows:

.. math::
    :label: gibbs_cost_function

    \mathrm{A_{system}(T, V)} & = \mathrm{G - PV}\\
                & = \sum_{i = 1}^{\mathrm{N_s + N_g}}
                \mu_i x_i - \mathrm{PV}


Where G is the Gibbs function. The system pressure P is in
bars. Dividing Eq. :eq:`gibbs_cost_function` by RT, we get,

.. math::
    :label: Atilde

    \mathrm{\tilde{A}(T, V)} = \frac{\mathrm{A_{system}}}{\mathrm{RT}}
    = \underbrace{\sum_{i = 1}^{\mathrm{N_s + N_g}}\tilde{\mu}_i x_i -
    \bar{x}^g}_{\mathrm{F(x)}}

Where :math:`\bar{x}^g` is the total number of moles of gas in the system.
:math:`\tilde{\mu}_i` is the reduced (dimension less) chemical
potential and can be given in terms of their standard chemical potential
functions as,

.. math::
    :label: chem_pot

    \tilde{\mu}_i =
    \begin{cases}
        {\left(\tilde{\mu}^o\right)}_i^g + \ln\mathrm{\frac{RT}{V}}
        + \ln x_i^g & \text{for }i = 1, 2 \ldots, \mathrm{N_g} \\
        {\left(\tilde{\mu}^{o}\right)}_{i}^{c} &
         \text{for }i = 1,2\ldots,\mathrm{N_s}
    \end{cases}


The superscript **g** and **c** are used for gaseous and condensed phase
species. For example, :math:`x_{i}^{g}` is the number of moles for
:math:`i^{th}` chemical species in gas phase.
:math:`{\left(\tilde{\mu}^o\right)}_{i}^{g}` is chemical potential
of :math:`i^{th}` gaseous chemical species in gas phase at standard
conditions. :math:`\mathrm{N_s}` is the total number of the condensed species,
:math:`\mathrm{N_g}` is the total number of the species in the gaseous phase.

While the free energy Eq. :eq:`Atilde` is minimized to solve for
equilibrium, the species have to satisfy non-negativity
constraint :math:`x_i\ge 0` and the constraint for element
conservation given as,

.. math::
    :label: constraint_equation

    \begin{align}
    \underbrace{\sum_{i = 1}^{\mathrm{N_g}} a_{ij}^{g}
    x_{i}^{g} + \sum _{i = 1}^{\mathrm{N_s}} a_{ij}^{c}x_{i}^{c} =
    b_j}_{\mathrm{C_j(x)}}            &&j = 1, 2, 3\ldots, \mathrm{N_e}
    \end{align}

Where :math:`a_{ij} ^{g}` is the number of atoms of :math:`j^{th}` element in
specie :math:`i` in the gas phase and :math:`a_{ij}^{c}` is the number of atoms
of :math:`j^{th}` element in :math:`i^{th}` species with
condensed phase. :math:`b_j` is the total number of moles of element j, originally present in system mixture.
The above free energy functions are minimised with two methods, viz., (1)
Quadratic gradient descent method (2) sequential least square minimisation
(SLSQP).

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Quadratic gradient descent method:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To minimise the free energy (Eq. :eq:`Atilde`) of the system
containing :math:`\mathrm{(x_1, x_2,\ldots x_{N_g},
\ldots,x_{N_g + N_s})}` moles of :math:`\mathrm{N_g}\) + \(\mathrm{N_s}` species, with
constraints described in Eq. :eq:`constraint_equation`, the method of
Lagrange :math:`'` s underdetermined multipliers [Bertsekas2014]_ [Weisstein2020]_ is used. Here, the formulation
is given for constant temperature and constant volume problem. Formulation for
constant pressure and constant temperature can be found in literature [White1958]_ [McBride_and_Gordon1996a]_ [Eriksson1979]_ [Eriksson_and_Hack1990]_.
The Lagrangian function to be minimized can be written as,

.. math::
    :label: eq_to_mini

    \mathrm{L} = \mathrm{\tilde{A}} - \sum_{j=1}^{N_e}\pi_j C_j(x)

Where, :math:`\pi_j` are the undetermined multipliers. The partial
derivatives of L with respect to :math:`i^{th}` chemical species can be given
after regrouping in terms of the gas species part and condensed species
part as,

.. math::
    :label: gas_part

    \frac{\partial L}{\partial x_i} = \mathrm{\frac{\partial L^g}
    {\partial x_i} +\frac{\partial L^c}{\partial x_i}}

Where, the :math:`\mathrm{\frac{\partial L^g}{\partial x_i}}` and
:math:`\mathrm{\frac{\partial L^c}{\partial x_i}}` are given by,


.. math::
    :label: gas_langrage_multiplier

    \begin{aligned}
    \mathrm{\frac{\partial L^g}{\partial x_i}} = (c_i + \ln x_i^g) -
    \sum_{j=1}^{\mathrm{N_e}}a_{ij}^{g}\pi_j = 0 && i = 1,2\ldots,\mathrm{N_g}
    \end{aligned}


.. math::
    :label: solid_langrage_multiplier

    \begin{aligned}
    \mathrm{\frac{\partial L^c}{\partial x_i}}
     = {\left(\tilde{\mu}^o\right)}_i^c - \sum_{j=1}^{\mathrm{N_e}}a_{ij}^{c}
    \pi_j = 0  &&i= 1,2\ldots,\mathrm{N_s}
    \end{aligned}

Where,

.. math::
    :label: c_i

    c_i = {\left(\tilde{\mu}^o\right)}_i^g + \ln\mathrm{\frac{RT}{V}}

To linearize Eq. :eq:`gas_langrage_multiplier`, Taylor series
expansion about :math:`y_{i}^{g}` is carried out, which results,

.. math::
    :label: taylor_expansion

    \begin{align}
    f_i + \frac{x_i^g}{y_i^g} - 1 - \sum_{j=1}^{\mathrm{N_e}}a_{ij}^{g}\pi_j
    = 0             &&i = 1,2\ldots,\mathrm{N_g}
    \end{align}

Where, :math:`f_{i}` can be given as,

.. math::
    :label: fi

    f_i = c_i + \ln y_i^g

From the Taylor expanded form, the improved mole numbers for the next
iteration can be obtained from the rearranged Eq. :eq:`taylor_expansion`

.. math::
    :label: improved_x

    x_i^g = -f_i y_i^g + y_i^g
    \left(\sum_{i=1}^{\mathrm{N_g}} \pi _{j} a_{ij}^{g} + 1\right)

Substituting Eq. :eq:`improved_x` in Eq. :eq:`constraint_equation` we have,

.. math::
    :label: r_matrix_equation

    \begin{align}
    \sum_{k=1}^{\mathrm{N_e}} r_{jk} \pi_k
    + \sum_{i=1}^{\mathrm{N_s}} a_{ij}^{c}x_{i}^{c} & = b_j
    + \sum_{i=1}^{\mathrm{N_g}} a_{ij}^{g} f_i y_i^g
    - \sum_{i=1}^{\mathrm{N_g}}a_{ij}^{g}y_i^g      &&j= 1,2\ldots, \mathrm{N_e}
    \end{align}

Where,

.. math::
    :label: r_ij

    \begin{align}
    r_{jk} & = r_{kj} = \sum_{i=1}^{\mathrm{N_g}}
    (a_{ij}^{g}a_{ik}^{g})y_{i}^{g} \text{          }j = 1, 2\ldots, \mathrm{N_e}
    \end{align}

The Eqs. :eq:`solid_langrage_multiplier` and :eq:`r_matrix_equation` are
solved simultaneously to get :math:`\pi _j` and :math:`x_i^c`. Using
:math:`\pi _j`, updated values of :math:`x^g` are obtained
using Eq. :eq:`improved_x`. Here, it should be noted that the formulation of
condensed species is such that, the moles of each condensed phase species is
directly obtained from the solution of Eq. :eq:`solid_langrage_multiplier` and
:eq:`r_matrix_equation`, without applying any correction.
If all :math:`x^g` values are positive, they are considered as the guessed value for the next iteration. If not, then guessed values are corrected with the following Eq.\ as [Gunnar_Eriksson1971]_,

.. math::
    :label: y_i_improved

    y_{i, improved}^{p} = y^{p}_{i} + \lambda (x^{p}_{i} - y^{p}_{i})


Where, p is the phase of chemical species (gaseous
species **g** or condensed species **c**). The :math:`\lambda` is
the correction factor, and can be given according
to [Gunnar_Eriksson1971]_,

.. math::
    :label: lambda

    \lambda = 0.99 \lambda ' (1 - 0.5 \lambda ')

Where, :math:`\lambda '` is the value required for the next step to
remain positive as given by,

.. math::
    :label: lambda_dash

    \lambda ' = \min_{i}\left(\frac{y^{g,c}_{i} }{(y^{g,c}_{i} - x^{g,c}_{i})}
    \right)

The corrected values of :math:`y^g` are considered for the next
iteration. Since, in the above formulation, the set of condensed species is
not known beforehand, for the first iteration, only gaseous species
equilibrium is obtained. In subsequent iterations, species are added such that
they reduce overall system :math:`'` s free energy. This can be assured by [McBride_and_Gordon1996a]_ [Gunnar_Eriksson1971]_,


.. math::
    :label: test_for_condesed_sp

    {\left(\tilde{\mu}^o\right)}_i^c -
    \sum_{j=1}^{\mathrm{N_e}}\pi_{j}b_{j} \le 0

Sometimes with particular species, the combination can lead to a set of the
dependent Eqs.\ in the formulation. Set of dependent species are identified
with the help of reduced row echelon form. Subsequently, dependent species, as
well as the matching combination of dependent species present in species set,
are removed temporarily, and only species which decrease overall free energy
of the system is included in species set. For example, species set containing,
U (L), :math:`\mathrm{UO_{2}(L)}` and :math:`\mathrm{U_4O_9(III)}`,
leads to linear dependence and all three species are removed, and species
which lower free energy of the system are included in species set.

The above Helmholtz function minimization subjected to mole number
conservation constraint is implemented in python and the code is hereafter
referred to as **MINICHEM** (MINImisation of CHEMical potentials). Once the equilibrium species moles in various
phases are determined, the release fraction of
the :math:`\mathrm{j^{th}}` chemical species are determined as follows:

.. math::

    \begin{aligned}
    \text{Elemental release fraction = }\mathrm{{RF(e)}_j} =
    \frac{{\sum\limits}_{i} \text{Number of moles of
    \(j^{th}\) element in \(i^{th}\) gaseous species}}{\text{Total mole
    inventory of \(j^{th}\) element}}
    \end{aligned}

The isotopic release fraction is defined as,

.. math::

    \text{Isotopic release fraction }{RF(iso)}_k = {\mathrm{RF(e)}_j} \times
    \text{isotopic fraction of \(k^{th}\) isotope}

