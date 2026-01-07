Input & Output Files
====================

The optimizer uses standard chemical file formats where possible, supplemented by a simple JSON schema for flexibility.

Input Files
-----------

Monomer XYZ File
~~~~~~~~~~~~~~~~

The primary input is a standard XYZ file specifying the geometry of a single monomer.

*   **Format:** Standard XYZ (number of atoms on first line, comment on second).
*   **Units:** Angstroms.
*   **Requirements:**
    *   Must be a valid geometry (no overlapping atoms).
    *   Orientation relative to axes does *not* matter; the program automatically aligns the molecular plane to the XY plane.

**Example (benzene.xyz):**

.. code-block:: text

   12
   Benzene molecule
   C         -1.21310        0.68840        0.00000
   C         -1.20280       -0.70630        0.00000
   C          0.01030       -1.39470        0.00000
   C          1.21310       -0.68840        0.00000
   C          1.20280        0.70630        0.00000
   C         -0.01030        1.39470        0.00000
   H         -2.15190        1.22980        0.00000
   H         -2.15130       -1.22940        0.00000
   H          0.01900       -2.48390        0.00000
   H          2.15190       -1.22980        0.00000
   H          2.15130        1.22940        0.00000
   H         -0.01900        2.48390        0.00000

Torsions Definition File (JSON)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For molecules with flexible bonds, you can specify which dihedrals to optimize using a JSON file.

**Schema:**

.. code-block:: json

   {
     "indexing": "0-based",
     "torsions": [
       {
         "name": "Linker Rotation",
         "atoms": [0, 1, 2, 3],
         "rotate_side": "d"
       }
     ]
   }

**Fields:**

*   ``indexing``: (Optional) ``"0-based"`` (default, Python/C style) or ``"1-based"`` (standard chemistry style).
*   ``torsions``: List of torsion definitions.
    *   ``atoms``: List of 4 atom indices defining the dihedral angle :math:`i-j-k-l`. The bond :math:`j-k` is the axis of rotation.
    *   ``rotate_side``: (Optional) Which side of the bond rotates.

        *   ``"d"`` (default): Distal side (finding atoms connected to atom :math:`l`).
        *   ``"p"``: Proximal side (finding atoms connected to atom :math:`i`).

    *   ``name``: (Optional) Label for logging.

Output Files
------------

optimization_results.txt
~~~~~~~~~~~~~~~~~~~~~~~~

The main summary file generated after a run. It contains:

1.  **Best Objective Value:** The optimized score (Energy + Penalties).
2.  **Binding Energy:** The pure binding energy (:math:`E_{bind}`) in kJ/mol.
3.  **Best Parameters:** The specific translation, rotation, and torsion angles that produced the best geometry.
4.  **Full Torsion Angles:** If torsions were used, the final values of all dihedrals are listed.

molecular_stack_Nmolecules.xyz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The final optimized geometry of the stack.

*   Contains :math:`N` copies of the monomer (e.g., 2 for a dimer).
*   This file is ready for visualization in VMD, PyMOL, or Mercury.

molecule_after_torsion.xyz
~~~~~~~~~~~~~~~~~~~~~~~~~~

The geometry of a **single monomer** in its optimized conformation (i.e., with the optimal internal torsion angles applied), before it is stacked. This is useful for inspecting how the monomer's shape changed during optimization.

traj.csv (Optional)
~~~~~~~~~~~~~~~~~~~

Generated if ``--print-trajectories`` is enabled. A raw CSV file containing the history of every evaluation.

**Columns:**
``iteration``, ``particle_id``, ``is_global_best``, ``c_theta``, ``s_theta``, ``Tx``, ``Ty``, ``Tz``, ``Cx``, ``Cy``, ``[torsions...]``, ``energy_score``
