# Molecular Stacking Optimization with PSO

This project implements a **Particle Swarm Optimization (PSO)** algorithm to find the optimal piling geometry of molecular dimers (or n-mers) by minimizing their binding energy. It leverages **xTB (Extended Tight-Binding)** for fast generic quantum mechanical energy evaluation.

## Features

- **Particle Swarm Optimization**: Custom implementation of PSO to explore the conformational space of molecular dimers.
- **xTB Integration**: Uses `xtb` (GFN2-xTB by default) to calculate single-point energies for objective function evaluation.
- **Geometry Transformation**: parameterized geometry employing 7 degrees of freedom:
    - Rotation: `cos_like`, `sin_like` (in-plane rotation)
    - Translation: `Tx`, `Ty`, `Tz`
    - Pivot Point: `Cx`, `Cy`
- **Multi-threading**: Supports parallel evaluation of the swarm using `ThreadPoolExecutor`.
- **Soft Clash Penalty**: Includes a smooth steric penalty to dampen high-energy clashes and improve stabilizer convergence.

## Requirements

### External Software
- **xTB**: The `xtb` executable must be installed and available in your system path.
    - [xTB GitHub Repository](https://github.com/grimme-lab/xtb)
    - [xTB Documentation](https://xtb-docs.readthedocs.io/en/latest/)

### Python Dependencies
Install the required packages using pip:

```bash
pip install -r requirements.txt
```

*Key dependencies:*
- `numpy`
- `ase` (Atomic Simulation Environment)

## Usage

1.  **Prepare Input**: Ensure you have an XYZ file for your monomer (e.g., `alkyl-pdi.xyz`).
2.  **Configure**: Edit `input.json` or modify the default config in `main.py` to set PSO parameters, swarm size, and file paths.
3.  **Run**:

```bash
python main.py
```

## Configuration

The simulation is controlled via `input.json`. Key sections include:

-   **`input_xyz`**: Path to the monomer XYZ file.
-   **`pso`**:
    -   `swarm_size`: Number of particles.
    -   `max_iters`: Maximum PSO iterations.
    -   `inertia`, `cognitive`, `social`: Hyperparameters for swarm dynamics.
-   **`tuning`**: options to enable automatic hyperparameter tuning before the main run.
-   **`xtb`**:
    -   `method`: xTB method (e.g., "gfn2", "gfnff").
    -   `charge`: Total system charge.
    -   `solvent`: Optional implicit solvent (e.g., "water").

## Outputs

-   **`outputs/`**: Directory containing optimized structures (e.g., `final_stack.xyz`).
-   **`energy_log.csv`**: Trace of the best global energy found at each iteration.
-   **`swarm_trace.csv`**: Detailed log of swarm positions and energies (if enabled).
-   **`optimization_results.txt`**: Summary of the best parameters and binding energy.
 Molecular Stacking Optimization

This project implements a Particle Swarm Optimization (PSO) algorithm to find the optimal stacking configuration of molecules using xTB (Extended Tight Binding) for energy evaluation.

## Features

- **Particle Swarm Optimization**: Optimizes 7 parameters (transformation matrix) to minimize binding energy.
- **xTB Integration**: Uses GFN2-xTB (via ASE or CLI) for accurate quantum mechanical energy calculations.
- **Parametric Transformation**: Defines molecular geometry using translation (Tx, Ty, Tz), rotation (Cx, Cy, θ), and pivot points.
- **N-Body Binding Energy**: Optimizes for stability in multi-layer stacks (not just dimers).
- **Hyperparameter Tuning**: Includes a random search tuner to find optimal PSO parameters (inertia, cognitive, social weights).

## Installation

1.  Clone the repository:
    ```bash
    git clone <your-repo-url>
    cd PSO
    ```

2.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

3.  Ensure `xtb` is installed and in your PATH (or use the ASE interface if configured).

## Usage

Run the main optimization script:

```bash
python main.py
```

Or with a custom configuration file:

```bash
python main.py my_config.json
```

## Documentation

Full documentation is available at: [Link to your ReadTheDocs here]

## Outputs

The script generates the following outputs in the `outputs/` directory:
- `energy_log.csv`: Iteration-level statistics.
- `swarm_trace.csv`: Detailed trajectory of every particle.
- `optimization_results.txt`: Summary of the best solution and parameters.
- `molecular_stack_Nmolecules.xyz`: The final optimized geometry.
