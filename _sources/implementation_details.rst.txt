Implementation Details
======================

This section provides a deep dive into the internal codebase of the π-Stack Optimizer. It is intended for developers or advanced users who want to understand *how* the code works under the hood.

Codebase Architecture
---------------------

The project is structured as a modular Python application. The core logic handles the orchestration between global optimization algorithms, geometric manipulations, and quantum chemical energy evaluations.

.. mermaid::

   graph TD
       A[pi-stack-generator.py] -->|Configures| B(BatchObjective)
       A -->|Selects| C(Optimizer)
       C -->|Queries| B
       B -->|Manages| D[XTBWorkerPool]
       D -->|Spawns| E[Worker Processes]
       E -->|Executes| F[xTB Binary]

Entry Point (`pi-stack-generator.py`)
-------------------------------------

The `main()` function in ``pi-stack-generator.py`` serves as the conductor:

1.  **Initialization**:
    *   Parses command line arguments (`argparse`).
    *   Sets up logging (`modules/logging_helpers.py`).
2.  **Preprocessing**:
    *   Reads the monomer XYZ file.
    *   **Validation**: Calls ``report_geometry_validation`` to ensure the input monomer has sane bond lengths and no overlaps.
    *   **Alignment**: Uses ``align_pi_core_to_xy`` to rotate the flat part of the molecule to the XY plane. This makes the optimization variables (translation/rotation) more physically intuitive.
    *   **Torsions**: Loads torsion definitions and optionally detects symmetry (`modules/torsion.py`).
3.  **Setup**:
    *   Initializes the ``BatchObjective`` class, which encapsulates the cost function.
    *   Initializes the ``XTBWorkerPool`` to manage parallel XTB processes.
    *   Selects the optimization strategy (PSO, GA, GWO) via ``create_optimizer``.
4.  **Execution List**:
    *   Runs the optimization loop: ``optimizer_obj.optimize(objective)``.
5.  **Result Handling**:
    *   Computes the final binding energy of the best solution.
    *   Writes output XYZ files and ``optimization_results.txt``.

Core Components
---------------

BatchObjective (`modules/objective.py`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This class defines the "landscape" that the optimizer explores. Its primary method is ``batch_evaluate(P)``, which takes a batch of parameter vectors and returns their energy scores.

**The Evaluation Pipeline:**

1.  **Parameter Decoding**:
    *   The input vector `P` contains rigid body parameters (rotation, translation) and torsion angles.
    *   Torsions are mapped from reduced space (if symmetry is used) to full space.

2.  **Geometry Reconstruction**:
    *   **Torsions**: ``_apply_torsions`` rotates the monomer fragments based on the angles.
    *   **Topology Check**: It validates that rotating bonds didn't accidentally break or form new bonds (e.g., steric clashes causing ring spearing).
    *   **Stack Building**: ``build_n_layer_stack`` creates the dimer (or n-mer) by applying the rigid body transformation to the monomer copies.

3.  **Fast Pre-Screening (The Penalties)**:
    *   Before running expensive quantum checks, the code runs fast geometric checks.
    *   **Clash Penalty**: Checks if atoms between layers are too close (< 1.6 Å).
    *   **Intramolecular Penalty**: Checks if the monomer itself is crashing into itself.
    *   *Optimization Note:* If a severe clash is found, the xTB calculation is **skipped entirely**, and a high penalty score is returned. This saves massive amounts of compute time.

4.  **Energy Calculation**:
    *   If the geometry is safe, a task is submitted to the ``XTBWorkerPool``.
    *   The binding energy is calculated as: :math:`E_{bind} = E_{stack} - N \times E_{monomer}`.

Geometry Engine (`modules/geometry.py`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This module handles the 3D math.

*   **Principal Component Analysis (PCA)**: Used in ``align_pi_core_to_xy``. It finds the plane of "flatness" in a molecule by analyzing the covariance matrix of atomic coordinates.
*   **Collision Detection**:
    *   Uses ``scipy.spatial.cKDTree`` (if available) for O(N log N) neighbor search.
    *   Falls back to NumPy broadcasting (O(N^2)) for environments without SciPy.

Concurrency Model (`modules/xtb_workers.py`)
--------------------------------------------

The software uses a **Process-based Parallelism** model via Python's `multiprocessing` module.

*   **Why Processes?** Python's Global Interpreter Lock (GIL) prevents threads from utilizing multiple CPU cores for CPU-bound tasks. Since we are launching external subprocesses (xTB), threads could work, but processes allow us to isolate the memory environment of each worker.
*   **The Pool**:
    *   The main process creates a queue of tasks.
    *   $N$ worker processes consume tasks from the queue.
    *   Each worker creates a unique temporary directory (`/tmp/xtb_run_...`) to run xTB without file collisions.
    *   Results are sent back via a results queue.

Optimization Algorithms (`modules/optimizer.py`)
------------------------------------------------

The optimizers are implemented as classes that share a common interface.

*   **PSO (Particle Swarm)**: Particles remember their personal best (`pbest`) and know the global best (`gbest`). Velocity updates mix these influences.
*   **GA (Genetic Algorithm)**: Uses tournament selection (pick random subset, choosing best) and simulated binary crossover (SBX) to breed new solutions.
*   **Hybrid PSO-NM**: Runs PSO for exploration, then periodically freezes the swarm to run a local Nelder-Mead simplex optimization on the best particle to "refine" the minimum.

Directory Structure Reference
-----------------------------

.. code-block:: text

   output/
   ├── modules/
   │   ├── geometry.py       # 3D math, PCA alignment, collision checks
   │   ├── objective.py      # The cost function (Geometry -> Energy)
   │   ├── optimizer.py      # PSO, GA, GWO algorithms
   │   ├── torsion.py        # Dihedral angle manipulation
   │   └── xtb_workers.py    # Parallel execution of xTB binaries
   ├── pi-stack-generator.py # Main entry point (Argument parsing, Orchestration)
   └── docs/                 # Documentation (Sphinx)
