Theoretical Background
======================

The **π-Stack Optimizer** solves the problem of finding energetically favorable configurations for molecular stacks. This is achieved by combining quantum chemical energy evaluations with global optimization algorithms.

Parameter Space Definition
--------------------------

The geometry of the molecular stack is defined by a set of parameters that describe the rigid body transformation of the monomer and its internal conformational flexibility. The total parameter space has dimensions :math:`7 + N`, where :math:`N` is the number of flexible dihedral angles.

Rigid Body Parameters (7 Dimensions)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The relative position and orientation of the monomer unit in the stack are described by:

1.  **Rotation Angle (:math:`\theta`):** To avoid discontinuities and periodic boundary issues associated with angles (e.g., :math:`0 \approx 2\pi`), the rotation is encoded using its sine and cosine components:
    
    .. math::

       c_\theta = \cos(\theta), \quad s_\theta = \sin(\theta)

    The optimization algorithms evolve :math:`c_\theta` and :math:`s_\theta` independently, and the actual angle is reconstructed as :math:`\theta = \text{atan2}(s_\theta, c_\theta)`.

2.  **Translation Vector (:math:`\mathbf{T}`):** Three parameters describe the translation of the monomer along the Cartesian axes:
    
    .. math::

       \mathbf{T} = [T_x, T_y, T_z]

3.  **Rotation Center (:math:`\mathbf{C}`):** Two parameters define the effective center of rotation in the X-Y plane (perpendicular to the stacking axis, usually Z):

    .. math::

       \mathbf{C} = [C_x, C_y]

    Combined with a fixed :math:`C_z = 0`, this allows for off-axis stacking arrangements.

Flexible Parameters (:math:`N` Dimensions)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For molecules with internal flexibility, specific rotatable bonds can be defined. Each bond :math:`i` contributes one degree of freedom:

.. math::

   \boldsymbol{\tau} = [\tau_1, \tau_2, \dots, \tau_N]

These are dihedral angles that modify the monomer's internal structure before the rigid body transformation is applied.

Objective Function
------------------

The goal of the optimization is to minimize the **binding energy** per monomer in the stack.

Energy Calculation
~~~~~~~~~~~~~~~~~~

The binding energy (:math:`E_{bind}`) is calculated using the **Semi-Empirical Tight-Binding (xTB)** method (specifically GFN2-xTB by default). For a system of :math:`N` monomers (typically a dimer, :math:`N=2`) with total energy :math:`E_{stack}` and isolated monomer energy :math:`E_{mono}`:

.. math::

   E_{bind} = \frac{E_{stack} - N \cdot E_{mono}}{N - 1}

Penalty Functions
~~~~~~~~~~~~~~~~~

To ensure physically realistic and computationally viable configurations, several penalty terms are added to the energy:

1.  **Hard Clash Penalty:** If any two atoms between adjacent monomers are closer than a cutoff distance (:math:`d_{inter} < 1.6 \text{ \AA}`):

    .. math::

       P_{clash} = 1000.0

2.  **Intramolecular Clash Penalty:** For flexible molecules, if atoms within the monomer overlap (:math:`d_{intra} < 1.2 \text{ \AA}`):
    
    .. math::

       P_{intra} = 1000.0

3.  **SCF Convergence Failure:** If the xTB calculation fails to converge (e.g., due to extreme geometry distortion):

    .. math::

       P_{fail} = 1,000,000

The final objective function value is:

.. math::

   f(\mathbf{x}) = E_{bind} + P_{clash} + P_{intra} + P_{fail}
