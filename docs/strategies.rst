Optimization Strategies
=======================

The π-Stack Optimizer implements several nature-inspired global optimization algorithms. This section details their mechanics and offers practical advice for tuning them.

Particle Swarm Optimization (PSO)
---------------------------------

PSO simulates a flock of birds (swarm) moving through the parameter space. It is the default algorithm because it offers a good balance between exploration (searching new areas) and exploitation (refining known good areas).

**Update Equations**

.. math::

   v_i^{t+1} = \omega v_i^t + c_1 r_1 (pbest_i - x_i^t) + c_2 r_2 (gbest - x_i^t)

.. math::

   x_i^{t+1} = x_i^t + v_i^{t+1}

**Parameter Tuning Guide**

*   **Inertia** (:math:`\omega`): Controls momentum.
    *   *High (e.g., 0.9)*: Particles fly through the space quickly. Good for initial exploration or large, flat energy landscapes.
    *   *Low (e.g., 0.4)*: Particles slow down and settle. Good for fine-tuning.
    *   *Default (0.73)*: A mathematically derived "magic number" that generally works well.

*   **Cognitive** (:math:`c_1`) vs **Social** (:math:`c_2`):
    *   *High Cognitive*: Particles are nostalgic; they return to *their own* best discoveries. Maintains diversity.
    *   *High Social*: Particles are conformist; they rush to the *swarm's* best discovery. Converges fast but risks premature convergence to local minima.
    *   *Recommendation*: Keep them equal (1.5) or slightly bias Social if you need faster results.

Genetic Algorithm (GA)
----------------------

The GA mimics biological evolution. It is often superior for "rough" energy landscapes with many sharp peaks, as it can jump across barriers via crossover.

**Mechanics**
1.  **Selection**: Uses Tournament Selection. A small group (size 3) is picked at random, and the fittest one becomes a parent.
2.  **Crossover**: Uses Simulated Binary Crossover (SBX). Parents mix their parameters to create offspring.
3.  **Mutation**: Adds Gaussian noise to parameters. This prevents the population from becoming identical.

**Tuning Tips**
*   **Population Size**: Larger is better for complex molecules (flexible linkers). Use 100+ for difficult cases.
*   **Mutation Rate**: If the algorithm gets stuck early, increase this (e.g., to 0.2).

Grey Wolf Optimizer (GWO)
-------------------------

GWO simulates the leadership hierarchy of wolves. It is mathematically unique because it tracks the **top 3** solutions (:math:`\alpha, \beta, \delta`) simultaneously, rather than just one global best.

**When to use GWO?**
*   GWO is often **faster** than PSO for simple rigid-body stacking (dimers without torsions).
*   It has fewer parameters to tune (no inertia or coefficients), making it easier to use "out of the box."

Hybrid PSO + Nelder-Mead
------------------------

This is the most expensive but most accurate method.

1.  **Stage 1 (Global)**: PSO explores the landscape to find a promising basin of attraction.
2.  **Stage 2 (Local)**: Periodically, the Nelder-Mead simplex algorithm is launched from the best particle to "drill down" to the exact bottom of the well.

**Best For:**
*   Final production runs where you need the energy accurate to 0.01 kJ/mol.
*   Situations where the energy landscape is "pitted" (many small local minima).


