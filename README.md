# PSO Molecular Stacking Optimization

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
# PSO
