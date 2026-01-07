User Guide
==========

This guide provides a comprehensive reference for using the **π-Stack Optimizer**. The core interface is the command-line tool ``pi-stack-generator``, which orchestrates the entire workflow from geometry preprocessing to global optimization and local refinement.

Command Line Interface (CLI)
----------------------------

The general syntax is:

.. code-block:: bash

   pi-stack-generator [xyz_file] [options]


Positional Arguments
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Argument
     - Description
   * - ``xyz_file``
     - Path to the input monomer structure in XYZ format. The geometry should be reasonable (no severe clashes), but orientation does not matter as the program will automatically alignment it to the XY plane.

General Options
~~~~~~~~~~~~~~~

.. option:: --n-layer <int>

   Number of molecular layers in the stack. Default is ``2`` (dimer).
   
   *Example:* ``--n-layer 3`` for a trimer optimization.

.. option:: --optimizer <str>

   The global optimization algorithm to use.
   
   *   ``pso``: Particle Swarm Optimization (Default). Robust and general-purpose.
   *   ``ga``: Genetic Algorithm. Good for discrete-like landscapes.
   *   ``gwo``: Grey Wolf Optimizer. Often converges faster but can be less exhaustive.
   *   ``pso-nm``: Hybrid PSO + Nelder-Mead. Best for high-precision refinement.

.. option:: --max-iters <int>

   Maximum number of iterations (or generations for GA). Default is ``300``.

.. option:: --seed <int>

   Random seed for reproducibility. Default is ``42``.

Optimization Parameters
~~~~~~~~~~~~~~~~~~~~~~~

**Particle Swarm (PSO)**

.. option:: --swarm-size <int>
   
   Number of particles in the swarm. Default: ``60``.

.. option:: --inertia <float>
   
   Inertia weight (:math:`\omega`). Controls how much the particle keeps its previous velocity. Higher values (>0.7) facilitate exploration; lower values (<0.5) facilitate exploitation. Default: ``0.73``.

.. option:: --cognitive <float>
   
   Cognitive coefficient (:math:`c_1`). Weight for the particle's *own* best known position. Default: ``1.50``.

.. option:: --social <float>
   
   Social coefficient (:math:`c_2`). Weight for the *global* best known position. Default: ``1.50``.

**Genetic Algorithm (GA)**

.. option:: --ga-population <int>
   
   Size of the population. Default: ``80``.

.. option:: --ga-mutation-rate <float>
   
   Probability of mutation per parameter gene. Default: ``0.10``.

xTB Backend Options
~~~~~~~~~~~~~~~~~~~

.. option:: --workers <int>

   Number of parallel worker processes to spawn. Set this to match your CPU core count. Default: ``4``.

.. option:: --method <str>

   xTB method to use. Options: ``gfn2`` (default, best accuracy), ``gfn1``, ``gfn0`` (fastest), ``gfnff``.

.. option:: --solvent <str>

   Implicit solvent model (e.g., ``benzene``, ``water``). Default: None (gas phase).

Usage Examples
--------------

Basic Dimer Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

The simplest usage requires only an XYZ file. This will perform a dimer optimization using PSO with default settings.

.. code-block:: bash

   pi-stack-generator monomer.xyz

High-Accuracy Trimer Search
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To optimize a trimer (3-stack) using the hybrid PSO-NM algorithm for higher precision, and increasing the swarm size for better coverage:

.. code-block:: bash

   pi-stack-generator monomer.xyz --n-layer 3 --optimizer pso-nm --swarm-size 100 --max-iters 500

Speed-Optimized Screening
~~~~~~~~~~~~~~~~~~~~~~~~~

For quick screening of many molecules, use the GFN0-xTB method (faster but less accurate) and the GWO optimizer (fast convergence):

.. code-block:: bash

   pi-stack-generator monomer.xyz --method gfn0 --optimizer gwo --max-iters 100

Handling Flexible Molecules
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your molecule has rotatable bonds (dihedrals) that should be part of the optimization, define them in a JSON file and pass it via ``--torsions-file``.

.. code-block:: bash

   pi-stack-generator monomer.xyz --torsions-file torsions.json --enable-symmetric-torsions

See :doc:`inputs_outputs` for the format of the torsions file.

Restarting / Logging
--------------------

The optimizer prints progress to specific log files:

*   ``optimization_results.txt``: Final summary.
*   ``traj.csv``: (Optional) Full trajectory of the optimization if ``--print-trajectories`` is used.

To track progress in real-time:

.. code-block:: bash

   tail -f optimization_results.txt
