Usage
=====

The π-Stack Optimizer provides a command-line interface for generating molecular stacks and performing global optimization.

Basic Usage
-----------

The primary command for running optimizations is ``pi-stack-generator`` (if using the wrapper) or by executing the python module directly.

.. code-block:: bash

   pi-stack-generator [INPUT_FILE] [OPTIONS]

Arguments
---------

.. option:: input_file

   Path to the input file (JSON or XYZ format) defining the molecular system.

.. option:: --optimizer, -opt

   The optimization algorithm to use. Available options:
   
   * ``pso``: Particle Swarm Optimization (Default)
   * ``ga``: Genetic Algorithm
   * ``gwo``: Grey Wolf Optimizer
   * ``pso_nm``: Hybrid PSO + Nelder-Mead

.. option:: --workers, -w

   Number of parallel workers for energy evaluation.

.. option:: --config, -c

   Path to a conceptual configuration file for advanced parameter tuning.

Examples
--------

**Run a standard PSO optimization:**

.. code-block:: bash

   pi-stack-generator benzene_dimer.xyz --optimizer pso

**Run with Genetic Algorithm and 4 workers:**

.. code-block:: bash

   pi-stack-generator input.json --optimizer ga --workers 4

Hyperparameter Optimization
---------------------------

For advanced users, the ``pi-hyperopt`` tool allows for tuning the meta-parameters of the optimization algorithms using Optuna.

.. code-block:: bash

   pi-hyperopt --study-name mysstudy --n-trials 100
