Input and Output Files
======================

Input Files
-----------

XYZ File
~~~~~~~~

The primary input is a standard **XYZ file** containing the atomic coordinates of a single monomer unit.
*   **Format:**

    *   Line 1: Number of atoms.
    *   Line 2: Comment line (can include charge/multiplicity info).
    *   Line 3+: Element symbol and X, Y, Z coordinates.

**Example (benzene.xyz):**

.. code-block:: text

   12
   Benzene molecule
   C        0.00000        1.39700        0.00000
   C        1.21000        0.69800        0.00000
   ...

torsions.json
~~~~~~~~~~~~~

For molecules with flexible bonds, a JSON file defines the definition of rotatable dihedrals.

*   **Structure:** A list of torsion objects.
*   **Fields:**

    *   `indices`: Array of 4 atom indices (0-based or 1-based, system detects automatically) defining the dihedral angle.
    *   `rotate_side`: Specifies which part of the molecule moves (`"left"`, `"right"`, or `"center"`).

**Example:**

.. code-block:: json

   [
       {
           "indices": [0, 1, 4, 5],
           "rotate_side": "right"
       }
   ]

Output Files
------------

Generated upon completion of a run:

1.  **optimized_monomer_[name].xyz**
    
    The geometry of the single monomer unit in its optimized internal conformation (if flexible bonds were optimized).

2.  **optimized_stack_[name].xyz**
    
    A resulting stack structure (typically 10 layers) generated using the optimal rigid body parameters found. This is useful for visual inspection.

3.  **optimization_results.txt**
    
    A summary text file containing:
    *   Final optimized parameters (:math:`T, C, \theta`).
    *   Final Binding Energy.
    *   Details of the xTB calculation.

4.  **pso_trajectory.csv** (Optional)
    
    If trajectory logging is enabled, this file contains the history of the optimization, including:
    *   Iteration number.
    *   Particle ID.
    *   Parameter values.
    *   Energy/Fitness values.
