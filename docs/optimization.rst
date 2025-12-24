Optimization Algorithms
=======================

The π-Stack Optimizer implements several nature-inspired global optimization algorithms to explore the potential energy surface.

Particle Swarm Optimization (PSO)
---------------------------------

PSO simulates the behavior of a flock of birds. A swarm of particles moves through the parameter space, with each particle's velocity influenced by its own best known position and the swarm's global best position.

**Update Equations:**

For each particle $i$ at iteration $t$:

.. math::

   v_i^{t+1} = \omega v_i^t + c_1 r_1 (pbest_i - x_i^t) + c_2 r_2 (gbest - x_i^t)

.. math::

   x_i^{t+1} = x_i^t + v_i^{t+1}

Where:
*   $\omega$: Inertia weight (controls momentum).
*   $c_1, c_2$: Cognitive and Social coefficients.
*   $r_1, r_2$: Random numbers in $[0, 1]$.
*   $pbest_i$: Personal best position of particle $i$.
*   $gbest$: Global best position of the swarm.

Genetic Algorithm (GA)
----------------------

The GA evolves a population of solutions using operators inspired by biological evolution.

**Key Components:**
*   **Selection:** Tournament selection is used to choose parents for the next generation.
*   **Crossover:** Simulated Binary Crossover (SBX) creates offspring by recombining parent parameters.
*   **Mutation:** Gaussian mutation adds random perturbations to maintain diversity.
*   **Elitism:** A fraction of the best individuals are preserved unchanged.

Grey Wolf Optimizer (GWO)
-------------------------

GWO mimics the leadership hierarchy and hunting mechanism of grey wolves.

**Hierarchy:**
*   $\alpha$ (Alpha): The best solution found so far.
*   $\beta$ (Beta): The second best solution.
*   $\delta$ (Delta): The third best solution.
*   $\omega$ (Omega): The rest of the search agents.

**Position Update:**
The position of a wolf is updated by averaging the vectors towards the three leaders:

.. math::

   \vec{D}_\alpha = |\vec{C}_1 \cdot \vec{X}_\alpha - \vec{X}|, \quad \vec{X}_1 = \vec{X}_\alpha - \vec{A}_1 \cdot \vec{D}_\alpha

.. math::

   \vec{X}(t+1) = \frac{\vec{X}_1 + \vec{X}_2 + \vec{X}_3}{3}

PSO + Nelder-Mead Hybrid
------------------------

This hybrid approach combines the global search capability of PSO with the local refinement of the Nelder-Mead simplex method. 
Periodically, the global best solution found by PSO is used as a starting point for a local Nelder-Mead optimization run to fine-tune the minimum.

Hyperparameter Optimization
---------------------------

The optimizer includes a built-in hyperparameter tuning framework using **Optuna**. This allows for the automatic discovery of optimal algorithm parameters (e.g., PSO inertia, GA mutation rate) for a specific molecular system.

**Optimization Modes:**
*   **Sequential:** Optimizes parameters for one molecule at a time, using previous results to warm-start the search.
*   **Joint:** Finds a single robust set of hyperparameters that performs well across a dataset of multiple molecules.
