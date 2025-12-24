Usage and CLI Reference
=======================

The **π-Stack Optimizer** provides a robust command-line interface.

Main Command: pi-stack-generator
--------------------------------

Run the optimization using the wrapper script or python module.

.. code-block:: bash

   pi-stack-generator [INPUT_FILE] [OPTIONS]

General Options
~~~~~~~~~~~~~~~

.. option:: input_file

   (Required) Path to the `.xyz` or `.json` file describing the system.

.. option:: --optimizer, -opt

   Optimization algorithm to use.
   *   Default: ``pso``
   *   Choices: ``pso``, ``ga``, ``gwo``, ``pso_nm``

.. option:: --workers, -w

   Number of parallel xTB workers.
   *   Default: All available CPUs

.. option:: --max-iters, -n

   Maximum number of optimization iterations.
   *   Default: ``100``

.. option:: --seed

   Random seed for reproducibility.
   *   Default: ``None``

Algorithm Specific Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**PSO Options:**

.. option:: --swarm-size

   Number of particles in the swarm.
   *   Default: ``30``

.. option:: --inertia

   Inertia weight (:math:`\omega`).
   *   Default: ``0.73``

.. option:: --cognitive

   Cognitive coefficient (:math:`c_1`).
   *   Default: ``1.50``

.. option:: --social

   Social coefficient (:math:`c_2`).
   *   Default: ``1.50``

**Genetic Algorithm Options:**

.. option:: --ga-population

   Population size.
   *   Default: ``80``

.. option:: --ga-mutation-rate

   Probability of mutation.
   *   Default: ``0.10``

.. option:: --ga-crossover-rate

   Probability of crossover.
   *   Default: ``0.90``

**Grey Wolf Optimizer Options:**

.. option:: --gwo-pack-size

   Number of wolves (search agents).
   *   Default: ``50``

.. option:: --gwo-a-start

   Initial value of the convergence parameter :math:`a`.
   *   Default: ``2.0``

.. option:: --gwo-a-end

   Final value of the convergence parameter :math:`a`.
   *   Default: ``0.0``

Physical & backend Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: --r-min

   Minimum allowed intermolecular distance (Penalty threshold).
   *   Default: ``3.0`` Å

.. option:: --clash-cutoff

   Hard clash cutoff distance.
   *   Default: ``1.6`` Å

.. option:: --gfn

   xTB GFN parameterization.
   *   Default: ``2``

Hyperparameter Tuning (pi-hyperopt)
-----------------------------------

The ``pi-hyperopt`` command is used for automated parameter tuning using Optuna.

.. code-block:: bash

   pi-hyperopt --molecules-root ./molecules --trials-per-molecule 50 --joint

Arguments:

.. option:: --molecules-root

   Directory containing subfolders of molecules to train on.

.. option:: --trials-per-molecule

   Number of Optuna trials to run per molecule.
   *   Default: ``50``

.. option:: --joint

   If set, runs Joint Optimization mode to find a single best parameter set for the dataset.
