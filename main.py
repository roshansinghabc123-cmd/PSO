
from __future__ import annotations


import os
import re
import json
import time
import math
import shutil
import hashlib
import tempfile
import subprocess
from dataclasses import dataclass, asdict
from typing import Callable, Dict, Optional, Tuple, List, Any

from concurrent.futures import ThreadPoolExecutor

import numpy as np


# ============================================================================
#                               Global constants
# ============================================================================

# Parameter naming & indexing (fixed ordering)
PARAM_NAMES = ["cos_like", "sin_like", "Tx", "Ty", "Tz", "Cx", "Cy"]
IDX = {name: i for i, name in enumerate(PARAM_NAMES)}

# Physical constants / conversions
HARTREE_TO_KJMOL = 2625.499638  # kJ/mol per Eh
EV_TO_HARTREE = 1.0 / 27.211386245988


# ============================================================================
#                           Transformation utilities
# ============================================================================

def angle_from_cos_sin_like(c_like: float, s_like: float) -> float:
    """
    Map unconstrained 'cos_like'/'sin_like' to a valid in-plane angle via atan2.
    This avoids having to normalize or enforce cos^2 + sin^2 = 1 explicitly.
    Returns the angle (radians).
    """
    return float(np.arctan2(s_like, c_like))


def transformation_matrix(params: np.ndarray) -> np.ndarray:
    """
    Build the 4x4 homogeneous transformation matrix for the 7-parameter system.

    The parameters vector ``params`` must contain 7 elements in the following order:
    ``[cos_like, sin_like, Tx, Ty, Tz, Cx, Cy]``.

    - ``cos_like``, ``sin_like``: Unconstrained parameters defining the rotation angle.
    - ``Tx``, ``Ty``, ``Tz``: Translation vector components.
    - ``Cx``, ``Cy``: Coordinates of the pivot point (z=0 implied).

    The transformation represents a rotation about the z-axis (in-plane) around the
    pivot point :math:`C = (Cx, Cy, 0)`, followed by a translation :math:`T = (Tx, Ty, Tz)`.

    Args:
        params (np.ndarray): A 1D array of shape (7,) containing the transformation parameters.

    Returns:
        np.ndarray: A 4x4 homogeneous transformation matrix that acts on *column* vectors.
                    (Note: this codebase typically applies it to row vectors via ``M.T``).

    Raises:
        ValueError: If ``params`` does not have the shape (7,).
    """
    if params.shape != (7,):
        raise ValueError("Expected params shape (7,) in order: "
                         "[cos_like, sin_like, Tx, Ty, Tz, Cx, Cy].")
    c_like, s_like, Tx, Ty, Tz, Cx, Cy = (float(params[i]) for i in range(7))

    theta = angle_from_cos_sin_like(c_like, s_like)
    c = float(np.cos(theta))
    s = float(np.sin(theta))

    # Δ = T + C - R*C
    dx = Tx + Cx - (Cx * c - Cy * s)
    dy = Ty + Cy - (Cx * s + Cy * c)
    dz = Tz

    M = np.array([
        [c, -s, 0.0, dx],
        [s,  c, 0.0, dy],
        [0.0, 0.0, 1.0, dz],
        [0.0, 0.0, 0.0, 1.0],
    ], dtype=float)
    return M


def apply_transform(coords: np.ndarray, M: np.ndarray) -> np.ndarray:
    """
    Apply a homogeneous transform matrix ``M`` to an array of Cartesian coordinates.

    Args:
        coords (np.ndarray): An (N, 3) array of coordinates [x, y, z].
        M (np.ndarray): A 4x4 homogeneous transformation matrix.

    Returns:
        np.ndarray: An (N, 3) array of transformed coordinates.

    Raises:
        ValueError: If ``coords`` is not an (N, 3) array.
    """
    arr = np.asarray(coords, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError("coords must be an (N,3) array.")
    ones = np.ones((arr.shape[0], 1), dtype=float)
    homo = np.hstack([arr, ones])          # (N,4)
    out = homo @ M.T                       # row-wise application to column-form M
    return out[:, :3]


# Optional helpers for bounds (kept inert unless bounds provided in config)
def _apply_bounds_in_place(vec: np.ndarray,
                           bounds: Dict[str, Tuple[float, float]]) -> None:
    for name, (lo, hi) in bounds.items():
        i = IDX[name]
        vec[i] = np.clip(vec[i], lo, hi)


def _sample_within_bounds(rng: np.random.Generator,
                          n: int,
                          bounds: Dict[str, Tuple[float, float]]) -> np.ndarray:
    dim = len(PARAM_NAMES)
    X = np.empty((n, dim), dtype=float)
    for name, i in IDX.items():
        if name in bounds:
            lo, hi = bounds[name]
            X[:, i] = rng.uniform(lo, hi, size=n)
        else:
            X[:, i] = rng.normal(0.0, 1.0, size=n)
    return X


# ============================================================================
#                           Soft clash/penalty term
# ============================================================================

def clash_penalty(coords_a: np.ndarray,
                  coords_b: np.ndarray,
                  cutoff: float = 1.6) -> float:
    """
    Compute a smooth overlap penalty between two sets of coordinates (monomers A and B).

    This function serves as a stabilizer to discourage very small interatomic distances
    which can cause SCF convergence failures in xTB. It is not a physical repulsion model.

    Args:
        coords_a (np.ndarray): Coordinates of monomer A (Na, 3).
        coords_b (np.ndarray): Coordinates of monomer B (Nb, 3).
        cutoff (float, optional): Distance scaling factor for the exponential penalty. Defaults to 1.6.

    Returns:
        float: A scalar penalty value (unitless). Larger values indicate significant overlap.
    """
    A = np.asarray(coords_a, dtype=float)
    B = np.asarray(coords_b, dtype=float)
    diff = A[:, None, :] - B[None, :, :]          # (Na, Nb, 3)
    d = np.linalg.norm(diff, axis=-1) + 1e-12
    return float(np.sum(np.exp(-(d / cutoff) ** 2)))


# ============================================================================
#                               PSO structures
# ============================================================================

@dataclass
class PSOConfig:
    """
    Configuration parameters for the Particle Swarm Optimization algorithm.

    Attributes:
        swarm_size (int): Number of particles in the swarm.
        max_iters (int): Maximum number of iterations to run.
        inertia (float): Inertia weight (w), controlling the impact of previous velocity.
        cognitive (float): Cognitive coefficient (c1), pulling particles to their personal best.
        social (float): Social coefficient (c2), pulling particles to the global best.
        seed (Optional[int]): Random seed for reproducibility.
        verbose_every (int): Frequency (in iterations) of progress printing.
        early_stop_patience (int): Iterations to wait for improvement before stopping early.
        early_stop_min_delta (float): Minimum improvement required to reset patience.
        check_every (int): Frequency (in iterations) to check for early stopping.
        pso_workers (int): Number of parallel workers for objective evaluation.
    """
    # Swarm + dynamics
    swarm_size: int = 80
    max_iters: int = 500
    inertia: float = 0.7298
    cognitive: float = 1.49618
    social: float = 1.49618
    # Reproducibility
    seed: Optional[int] = 42
    # Verbosity / logging cadence (to stdout and energy_log.csv)
    verbose_every: int = 10
    # Early stopping
    # stop if no improvement for this many iterations
    early_stop_patience: int = 50
    # min improvement (kJ/mol) to reset patience
    early_stop_min_delta: float = 1e-2
    check_every: int = 5                  # evaluate the stopping rule every N iterations
    # Parallel evals of the objective (each *calls* xTB which itself can use threads)
    pso_workers: int = 1                  # keep 1 by default to avoid oversubscription


class PSO:
    """
    Particle Swarm Optimization (PSO) implementation.

    This class manages the swarm state, updates particle positions and velocities,
    and tracks global and personal bests.
    """

    def __init__(self, cfg: PSOConfig, rng: Optional[np.random.Generator] = None):
        """
        Initialize the PSO optimizer.

        Args:
            cfg (PSOConfig): Configuration settings for the PSO.
            rng (Optional[np.random.Generator]): Pre-seeded random number generator.
                                                 If None, a new one is created using ``cfg.seed``.
        """
        self.cfg = cfg
        self.dim = len(PARAM_NAMES)
        self.rng = rng if rng is not None else np.random.default_rng(cfg.seed)

    # ---------------------- Initialization ----------------------
    def _initialize_swarm(self,
                          bounds: Optional[Dict[str,
                                                Tuple[float, float]]] = None
                          ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Initialize particle positions P and velocities V.
        If `bounds` is provided, sample uniform in bounds; otherwise use
        a reasonable unconstrained Gaussian initialization.
        """
        if bounds:
            P = _sample_within_bounds(self.rng, self.cfg.swarm_size, bounds)
        else:
            P = np.zeros((self.cfg.swarm_size, self.dim), dtype=float)
            P[:, IDX["cos_like"]] = self.rng.normal(
                0.0, 1.0, size=self.cfg.swarm_size)
            P[:, IDX["sin_like"]] = self.rng.normal(
                0.0, 1.0, size=self.cfg.swarm_size)
            for name, sd in (("Tx", 2.0), ("Ty", 2.0), ("Cx", 2.0), ("Cy", 2.0)):
                P[:, IDX[name]] = self.rng.normal(
                    0.0, sd, size=self.cfg.swarm_size)
            # Positive z-separation to start around ~π-stack distances
            P[:, IDX["Tz"]] = self.rng.normal(
                3.5, 1.0, size=self.cfg.swarm_size)

        V = self.rng.normal(0.0, 0.5, size=P.shape)
        return P, V

    # ----------------------- Main optimize ----------------------
    def optimize(self,
                 score_fn: Callable[[np.ndarray], float],
                 bounds: Optional[Dict[str, Tuple[float, float]]] = None,
                 on_iter: Optional[Callable[[int, np.ndarray, np.ndarray, np.ndarray],
                                            None]] = None
                 ) -> Tuple[np.ndarray, float]:
        """
        Run the PSO and return (gbest_params, gbest_value).

        Args:
          score_fn: maps a single parameter vector → scalar loss (to be *minimized*).
          on_iter : optional callback executed each iteration with signature
                    on_iter(it, P, V, fitness). Use for logging/tracing.

        Returns:
          best_params, best_value
        """
        cfg = self.cfg
        P, V = self._initialize_swarm(bounds)

        # Helper to evaluate a batch (vectorized over particles) with optional threads.
        def evaluate_batch(batch: np.ndarray) -> np.ndarray:
            if cfg.pso_workers <= 1:
                return np.array([score_fn(p) for p in batch], dtype=float)
            else:
                with ThreadPoolExecutor(max_workers=cfg.pso_workers) as ex:
                    return np.array(list(ex.map(score_fn, list(batch))), dtype=float)

        # Evaluate initial fitness
        fitness = evaluate_batch(P)
        pbest = P.copy()
        pbest_val = fitness.copy()
        g_idx = int(np.argmin(pbest_val))
        gbest = pbest[g_idx].copy()
        gbest_val = float(pbest_val[g_idx])
        best_seen = gbest_val

        last_improve_it = 0  # iteration index when we last improved by min_delta

        # Main loop
        for it in range(1, cfg.max_iters + 1):
            r1 = self.rng.random(P.shape)
            r2 = self.rng.random(P.shape)

            # Velocity & position update
            V = (cfg.inertia * V
                 + cfg.cognitive * r1 * (pbest - P)
                 + cfg.social * r2 * (gbest - P))
            P = P + V

            if bounds:
                # Hard-clip inside bounds (optional). This also dampens runaway particles.
                for i_name, (lo, hi) in bounds.items():
                    i = IDX[i_name]
                    P[:, i] = np.clip(P[:, i], lo, hi)

            # Evaluate
            fitness = evaluate_batch(P)

            # Update personal bests
            improved = fitness < pbest_val
            pbest[improved] = P[improved]
            pbest_val[improved] = fitness[improved]

            # Update global best
            g_idx = int(np.argmin(pbest_val))
            if pbest_val[g_idx] < gbest_val:
                gbest_val = float(pbest_val[g_idx])
                gbest = pbest[g_idx].copy()

            # Trace / logging hook
            if on_iter is not None:
                on_iter(it, P, V, fitness)

            # Verbose progress
            if cfg.verbose_every and (it % cfg.verbose_every == 0 or it == 1):
                print(f"[PSO] iter={it:4d}  best={gbest_val:.6f}")

            # Early stopping check (every `check_every` iterations)
            if (it % cfg.check_every) == 0:
                if (best_seen - gbest_val) > cfg.early_stop_min_delta:
                    best_seen = gbest_val
                    last_improve_it = it
                elif (it - last_improve_it) >= cfg.early_stop_patience:
                    print(f"[PSO] Early stop at iter={it} (no ≥{cfg.early_stop_min_delta:.3g} "
                          f"improvement for {cfg.early_stop_patience} iters).")
                    break

        return gbest, gbest_val


# ============================================================================
#                                 XYZ I/O
# ============================================================================

def read_xyz_file(filename: str) -> Tuple[np.ndarray, List[str]]:
    """
    Read XYZ format file and return coordinates and atom symbols.
    Assumes strict standard XYZ with correct atom count in line 1.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    n_atoms = int(lines[0].strip())
    atoms: List[str] = []
    coords: List[List[float]] = []
    for i in range(2, 2 + n_atoms):  # Skip header lines
        parts = lines[i].strip().split()
        atoms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(coords, dtype=float), atoms


def write_xyz_file(filename: str, coords: np.ndarray, atoms: List[str], comment: str = "") -> None:
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for atom, coord in zip(atoms, coords):
            f.write(
                f"{atom:2s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}\n")


def build_xyz_string(atoms: List[str], coords: np.ndarray, comment: str = "") -> str:
    lines = [str(len(atoms)), comment]
    for atom, (x, y, z) in zip(atoms, coords):
        lines.append(f"{atom:2s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(lines) + "\n"


# ============================================================================
#                               xTB integration
# ============================================================================

@dataclass
class XTBConfig:
    charge: int = 0                      # total charge
    mult: int = 1                        # spin multiplicity (1 = closed shell)
    method: str = "gfn2"                 # 'gfn2'|'gfn1'|'gfn0'|'gfnff'
    solvent: Optional[str] = None        # e.g., "water" (ALPB)
    backend: str = "auto"                # 'auto'|'ase'|'cli'
    executable: str = "xtb"              # CLI executable name/path
    keep_workdir: bool = False           # keep temp files for inspection
    work_base: Optional[str] = None      # base directory for tmp folders

    def ase_method(self) -> str:
        mapping = {"gfn2": "GFN2-xTB", "gfn1": "GFN1-xTB",
                   "gfn0": "GFN0-xTB", "gfnff": "GFN-FF"}
        return mapping.get(self.method.lower(), "GFN2-xTB")


class XTBRunner:
    """
    Two backends:
      - ASE calculator ('ase.calculators.xtb.XTB') if installed.
      - CLI fallback calling `xtb coord.xyz --gfn N --sp --chrg Q --uhf U [--alpb solvent]`.

    The runner caches single-point energies by hashing the atomic identities and
    coordinates. This avoids redundant external calls when the same structure is
    evaluated repeatedly.
    """

    def __init__(self, cfg: XTBConfig):
        self.cfg = cfg
        self.backend = self._choose_backend(cfg.backend)
        self.cache: Dict[Tuple[str, int, int,
                               str, Optional[str], str], float] = {}

    # -------- Backend selection --------

    def _choose_backend(self, preference: str) -> str:
        if preference in ("ase", "auto"):
            try:
                import ase  # noqa: F401
                from ase.calculators.xtb import XTB as ASE_XTB  # noqa: F401
                self._ASE_AVAILABLE = True
                return "ase"
            except Exception:
                pass
        if preference in ("cli", "auto"):
            xtb_path = shutil.which(self.cfg.executable)
            if xtb_path is not None:
                self._XTB_PATH = xtb_path
                return "cli"
        raise RuntimeError(
            "No xTB backend available. Install ASE (`pip install ase`) or ensure the `xtb` CLI is in PATH."
        )

    # -------- Cache key --------

    @staticmethod
    def _coords_hash(atoms: List[str], coords: np.ndarray) -> str:
        m = hashlib.sha1()
        m.update((",".join(atoms)).encode())
        m.update(np.asarray(coords, dtype=np.float64).tobytes())
        return m.hexdigest()

    # -------- Public API --------

    def single_point_hartree(self, atoms: List[str], coords: np.ndarray) -> float:
        """
        Calculate the single-point total energy in Hartree (Eh).

        This method caches results based on atomic identities and coordinates to avoid
        redundant calculations.

        Args:
            atoms (List[str]): List of atom symbols (e.g., ["C", "H", ...]).
            coords (np.ndarray): (N, 3) array of Cartesian coordinates.

        Returns:
            float: Total energy in Hartree.
        """
        key = (self._coords_hash(atoms, coords),
               self.cfg.charge, self.cfg.mult, self.cfg.method, self.cfg.solvent, self.backend)
        if key in self.cache:
            return self.cache[key]

        if self.backend == "ase":
            e_h = self._ase_single_point(atoms, coords)
        else:
            e_h = self._cli_single_point(atoms, coords)

        self.cache[key] = e_h
        return e_h

    # -------- ASE backend --------

    def _ase_single_point(self, atoms: List[str], coords: np.ndarray) -> float:
        from ase import Atoms
        from ase.calculators.xtb import XTB as ASE_XTB

        uhf = max(self.cfg.mult - 1, 0)
        kwargs = dict(method=self.cfg.ase_method(),
                      charge=self.cfg.charge, uhf=uhf)
        if self.cfg.solvent:
            kwargs["solvent"] = self.cfg.solvent
        calc = ASE_XTB(**kwargs)

        a = Atoms(symbols=atoms, positions=coords)
        a.calc = calc
        e_ev = float(a.get_potential_energy())       # ASE returns eV
        return e_ev * EV_TO_HARTREE                  # convert to Eh

    # -------- CLI backend --------

    def _cli_single_point(self, atoms: List[str], coords: np.ndarray) -> float:
        def method_to_flag(m: str) -> str:
            m = m.lower()
            return {"gfn2": "2", "gfn1": "1", "gfn0": "0", "gfnff": "ff"}.get(m, "2")

        # Prepare work directory
        if self.cfg.keep_workdir:
            workdir = tempfile.mkdtemp(
                prefix="xtb_keep_", dir=self.cfg.work_base)
        else:
            workdir = tempfile.mkdtemp(prefix="xtb_", dir=self.cfg.work_base)

        coord_path = os.path.join(workdir, "coord.xyz")
        with open(coord_path, "w") as f:
            f.write(build_xyz_string(atoms, coords, comment="xTB single-point"))

        cmd = [
            getattr(self, "_XTB_PATH", "xtb"),
            coord_path,
            "--gfn", method_to_flag(self.cfg.method),
            "--sp",
            "--chrg", str(self.cfg.charge),
            "--uhf", str(max(self.cfg.mult - 1, 0)),
        ]
        if self.cfg.solvent:
            cmd += ["--alpb", self.cfg.solvent]

        # Run xTB (inherits environment, including OMP_NUM_THREADS set in main())
        run = subprocess.run(cmd, capture_output=True, text=True,
                             encoding='utf-8', errors='replace', cwd=workdir)
        stdout = run.stdout or ""
        stderr = run.stderr or ""
        out = stdout + "\n" + stderr

        # Parse energy in Eh
        m = re.search(r"(?i)total\s+energy\s*=\s*([-\d\.Ee+]+)\s*Eh", out)
        if not m:
            m = re.search(r"(?i)total\s+energy\s+([-\d\.Ee+]+)\s*Eh", out)
        if m:
            e_h = float(m.group(1))
        else:
            # Last scan for "E =  ... Eh"
            m = re.search(r"(?i)\bE\s*=\s*([-\d\.Ee+]+)\s*Eh", out)
            if m:
                e_h = float(m.group(1))
            else:
                if not self.cfg.keep_workdir:
                    shutil.rmtree(workdir, ignore_errors=True)
                raise RuntimeError(
                    "xTB energy not found in output. Last 2000 chars:\n" +
                    out[-2000:]
                )

        if not self.cfg.keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)

        return e_h


# ============================================================================
#                             Configuration handling
# ============================================================================

DEFAULT_CONFIG: Dict[str, Any] = {
    # --- General input/output ---
    "input_xyz": "alkyl-pdi.xyz",
    "output_dir": "outputs",
    # Build a full stack of this many molecules for the *final* export
    "build_n": 10,
    # Use this many molecules in the objective energy E_n (can be smaller
    # than build_n for speed). Must be >= 2.
    "objective_n": 2,

    # --- Threading / parallelisation ---
    # Threads for OpenMP/BLAS-backed libs (xTB, MKL, OpenBLAS, etc.)
    "threads": 1,
    # Concurrent objective evaluations per PSO iteration.
    # Keep at 1 unless you understand your CPU topology to avoid oversubscription.
    "pso_workers": 1,

    # --- xTB settings ---
    "xtb": {
        "charge": 0,
        "mult": 1,
        "method": "gfn2",
        "solvent": None,
        "backend": "auto",
        "executable": "xtb",
        "keep_workdir": False,
        "work_base": None
    },

    # --- PSO settings (used after tuning unless `tuning.enabled` is False) ---
    "pso": {
        "swarm_size": 80,
        "max_iters": 500,
        "inertia": 0.7298,
        "cognitive": 1.49618,
        "social": 1.49618,
        "seed": 42,
        "verbose_every": 10,
        "early_stop_patience": 50,
        "early_stop_min_delta": 1e-2,
        "check_every": 5,
    },

    # Optional bounds
    "bounds": None,

    # --- Tuning the PSO hyperparameters ---
    "tuning": {
        "enabled": True,
        # How many random configs to try
        "num_candidates": 8,
        # Iterations per candidate during tuning (kept modest for speed)
        "iters_per_candidate": 120,
        # Evaluate each candidate across N different random seeds and average
        "n_seeds": 2,
        "random_seed": 2025,
        # Ranges to sample
        "ranges": {
            "inertia": [0.5, 0.9],
            "cognitive": [0.8, 2.5],
            "social": [0.8, 2.5],
            "swarm_size": [50, 140]  # integer range
        }
    }
}


def _deep_update(base: Dict[str, Any], updates: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge `updates` into `base` (in-place) and return `base`."""
    for k, v in updates.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _deep_update(base[k], v)
        else:
            base[k] = v
    return base


def load_config(path: Optional[str]) -> Dict[str, Any]:
    """
    Load a JSON config file if present; otherwise return defaults.
    Missing keys fall back to DEFAULT_CONFIG values.
    """
    cfg = json.loads(json.dumps(DEFAULT_CONFIG))  # deep copy
    if path and os.path.isfile(path):
        with open(path, "r") as f:
            user_cfg = json.load(f)
        _deep_update(cfg, user_cfg)
    else:
        # Also try a default 'input.json' if present
        if os.path.isfile("input.json"):
            with open("input.json", "r") as f:
                user_cfg = json.load(f)
            _deep_update(cfg, user_cfg)
    return cfg


def ensure_output_dir(d: str) -> str:
    os.makedirs(d, exist_ok=True)
    return d


# ============================================================================
#                         Helpers for parameter reporting
# ============================================================================

def params_to_dict(vec: np.ndarray) -> Dict[str, float]:
    return {name: float(vec[IDX[name]]) for name in PARAM_NAMES}


def summarize_params(vec: np.ndarray) -> str:
    d = params_to_dict(vec)
    theta = angle_from_cos_sin_like(d["cos_like"], d["sin_like"])
    deg = np.degrees(theta)
    return (
        f"θ = {deg:.3f}°  "
        f"Tx = {d['Tx']:.4f} Å,  Ty = {d['Ty']:.4f} Å,  Tz = {d['Tz']:.4f} Å,  "
        f"Cx = {d['Cx']:.4f} Å,  Cy = {d['Cy']:.4f} Å"
    )


# ============================================================================
#                        Objective: n-molecule binding energy
# ============================================================================

def build_stack_from_transform(monomer_coords_centered: np.ndarray,
                               atom_types: List[str],
                               M: np.ndarray,
                               n_layers: int) -> Tuple[np.ndarray, List[str]]:
    """
    Construct a stack of `n_layers` by repeatedly applying transform `M`
    to the centered monomer coordinates.
    """
    all_coords: List[np.ndarray] = []
    all_atoms: List[str] = []
    current = monomer_coords_centered.copy()
    for _ in range(n_layers):
        all_coords.append(current)
        all_atoms.extend(atom_types)
        current = apply_transform(current, M)
    return np.vstack(all_coords), all_atoms


def n_body_binding_energy_kj(vec: np.ndarray,
                             xtb: XTBRunner,
                             monomer_coords_centered: np.ndarray,
                             atom_types: List[str],
                             E_mono_h: float,
                             n_layers: int,
                             penalty_weight: float = 2.0,
                             clash_cutoff: float = 1.6) -> float:
    """
    Compute the *n-body* binding energy in kJ/mol:
        ΔE_bind(n) = (E_n - n * E_1) / (n - 1)
    using the xTB single-point energies. Adds a small smooth steric penalty.

    Notes:
      - When n=2, this reduces to the usual dimer binding energy.
      - Coordinates are constructed by repeated application of the same transform M.
    """
    # 1) Construct transform from parameters
    M = transformation_matrix(vec)

    # 2) Build an n-layer stack from a *centered* monomer
    coords_n, atoms_n = build_stack_from_transform(
        monomer_coords_centered, atom_types, M, n_layers=n_layers
    )

    # 3) Smooth steric penalty between *adjacent* layers (helps optimizer stability)
    penalty = 0.0
    if n_layers >= 2:
        layer_size = monomer_coords_centered.shape[0]
        prev = monomer_coords_centered
        current = apply_transform(monomer_coords_centered, M)
        penalty += clash_penalty(prev, current, cutoff=clash_cutoff)

    # 4) xTB energies
    try:
        E_n_h = xtb.single_point_hartree(atoms_n, coords_n)
    except Exception:
        # Return a large penalty if xTB fails (extreme overlaps etc.)
        return 1.0e6 + 1.0e4 * penalty

    # 5) Binding energy in kJ/mol
    dE_kj = ((E_n_h - n_layers * E_mono_h) *
             HARTREE_TO_KJMOL) / max(n_layers - 1, 1)
    return float(dE_kj + penalty_weight * penalty)


# ============================================================================
#                        PSO hyperparameter tuning (random)
# ============================================================================

@dataclass
class TuningConfig:
    enabled: bool = True
    num_candidates: int = 8
    iters_per_candidate: int = 120
    n_seeds: int = 2
    random_seed: int = 2025
    ranges: Dict[str, List[float]] = None  # set from DEFAULT_CONFIG

    def __post_init__(self):
        if self.ranges is None:
            self.ranges = {
                "inertia": [0.5, 0.9],
                "cognitive": [0.8, 2.5],
                "social": [0.8, 2.5],
                "swarm_size": [50, 140]
            }


def sample_hyperparameters(rng: np.random.Generator, ranges: Dict[str, List[float]]) -> Dict[str, Any]:
    """
    Uniform random sampling over provided ranges. Swarm size is integer.
    """
    def samp(a, b):
        return float(a + (b - a) * rng.random())

    hp = {
        "inertia": samp(*ranges["inertia"]),
        "cognitive": samp(*ranges["cognitive"]),
        "social": samp(*ranges["social"]),
        "swarm_size": int(round(samp(*ranges["swarm_size"])))
    }
    # Common good practice: keep c1 + c2 in ~[2, 4]; we already sample within 0.8..2.5
    return hp


def tune_pso_hyperparameters(base_pso_cfg: Dict[str, Any],
                             tcfg: TuningConfig,
                             build_objective: Callable[[PSOConfig], Callable[[np.ndarray], float]],
                             bounds: Optional[Dict[str, Tuple[float, float]]],
                             out_dir: str) -> Dict[str, Any]:
    """
    Random-search tuner: tries `tcfg.num_candidates` PSO settings (each averaged
    across `tcfg.n_seeds` random seeds for robustness) and selects the PSO config
    that attains the lowest objective value within a modest iteration budget.

    We log a CSV summary (`tuning_log.csv`) and return the best hyperparameters.
    """
    log_path = os.path.join(out_dir, "tuning_log.csv")
    with open(log_path, "w") as f:
        f.write(
            "candidate,seed,swarm_size,inertia,cognitive,social,iters,best_value_kjmol\n")

    rng = np.random.default_rng(tcfg.random_seed)

    best_hp: Dict[str, Any] = {}
    best_score = math.inf

    for cidx in range(1, tcfg.num_candidates + 1):
        hp = sample_hyperparameters(rng, tcfg.ranges)
        # Merge into a temporary PSOConfig
        tmp_cfg = dict(base_pso_cfg)
        tmp_cfg.update(hp)
        tmp_cfg["max_iters"] = int(tcfg.iters_per_candidate)
        # Keep verbosity low during tuning
        tmp_cfg["verbose_every"] = max(0, int(tcfg.iters_per_candidate // 4))

        seed_scores: List[float] = []
        for s in range(1, tcfg.n_seeds + 1):
            tmp_cfg["seed"] = int(10_000 * cidx + s)  # make distinct
            pso = PSO(PSOConfig(**tmp_cfg))
            # Build objective function bound to this PSO instance
            score_fn = build_objective(pso.cfg)
            # Run with minimal logging

            def on_iter(_it, _P, _V, _fit):  # no-op during tuning
                pass
            gbest, gbest_val = pso.optimize(score_fn, bounds, on_iter)
            seed_scores.append(float(gbest_val))

            with open(log_path, "a") as f:
                f.write(f"{cidx},{s},{tmp_cfg['swarm_size']},{tmp_cfg['inertia']:.6f},"
                        f"{tmp_cfg['cognitive']:.6f},{tmp_cfg['social']:.6f},"
                        f"{tmp_cfg['max_iters']},{gbest_val:.6f}\n")

        avg_score = float(np.mean(seed_scores))
        if avg_score < best_score:
            best_score = avg_score
            best_hp = hp

    # Return the best hyperparameters merged onto base config
    final_cfg = dict(base_pso_cfg)
    final_cfg.update(best_hp)
    return final_cfg


# ============================================================================
#                             Logging utilities
# ============================================================================

class EnergyLogger:
    """
    Simple CSV logger for iteration-level statistics (mirrored to stdout by PSO).
    Each row: iter, best, mean, std
    """

    def __init__(self, path: str):
        self.path = path
        with open(self.path, "w") as f:
            f.write("iter,best_kjmol,mean_kjmol,std_kjmol\n")

    def log(self, it: int, fitness: np.ndarray):
        best = float(np.min(fitness))
        mean = float(np.mean(fitness))
        std = float(np.std(fitness))
        with open(self.path, "a") as f:
            f.write(f"{it},{best:.6f},{mean:.6f},{std:.6f}\n")


class SwarmTraceLogger:
    """
    CSV logger for *every particle* at *every iteration*.
    Columns:
      iter, particle, fitness, cos_like, sin_like, Tx, Ty, Tz, Cx, Cy
    """

    def __init__(self, path: str):
        self.path = path
        with open(self.path, "w") as f:
            header = ["iter", "particle", "fitness_kjmol"] + PARAM_NAMES
            f.write(",".join(header) + "\n")

    def log(self, it: int, P: np.ndarray, fitness: np.ndarray):
        with open(self.path, "a") as f:
            for pid, (vec, fit) in enumerate(zip(P, fitness)):
                row = [str(it), str(
                    pid), f"{float(fit):.6f}"] + [f"{float(x):.10f}" for x in vec]
                f.write(",".join(row) + "\n")


def write_optimization_results(path: str,
                               best_params: np.ndarray,
                               best_val: float,
                               pso_cfg: PSOConfig,
                               objective_n: int) -> None:
    """
    Human- and machine-readable results file that **starts** with the marker
    'OPTIMIZATION RESULTS', so other programs can scan for it easily.
    """
    d = params_to_dict(best_params)
    theta = angle_from_cos_sin_like(d["cos_like"], d["sin_like"])
    deg = np.degrees(theta)

    with open(path, "w") as f:
        f.write("OPTIMIZATION RESULTS\n")
        f.write("# Best transform parameters (suitable for generating stacks)\n")
        for k in PARAM_NAMES:
            f.write(f"{k}: {d[k]:.10f}\n")
        f.write(f"angle_degrees: {deg:.6f}\n")
        f.write(f"objective_n_layers: {int(objective_n)}\n")
        f.write(f"best_objective_kjmol: {best_val:.6f}\n")
        f.write("\n# Recommended PSO hyperparameters (after tuning)\n")
        f.write(json.dumps(asdict(pso_cfg), indent=2) + "\n")


# ============================================================================
#                            Stack statistics helpers
# ============================================================================

def layer_centroid_z(coords_layer: np.ndarray) -> float:
    """Mean z of a single layer's coordinates."""
    return float(np.mean(coords_layer[:, 2]))


def stack_statistics_from_centroids(layer_z: List[float]) -> Dict[str, float]:
    """
    Compute robust stack statistics from per-layer z-centroids:
      - height = z_last - z_first
      - separations = diffs between consecutive z
      - average separation = mean(separations)
      - stdev separation = std(separations)
    """
    if len(layer_z) < 2:
        return {"stack_height": 0.0, "avg_sep": 0.0, "std_sep": 0.0}

    layer_z_sorted = sorted(layer_z)
    height = layer_z_sorted[-1] - layer_z_sorted[0]
    seps = np.diff(np.array(layer_z_sorted))
    return {
        "stack_height": float(height),
        "avg_sep": float(np.mean(seps)),
        "std_sep": float(np.std(seps)),
    }


# ============================================================================
#                                    main()
# ============================================================================

def main(config_path: Optional[str] = None) -> None:
    # -------------------------------
    # Load + prepare configuration
    # -------------------------------
    cfg = load_config(config_path)

    threads = int(cfg.get("threads", 1))
    for var in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "ASE_CPUS"]:
        os.environ[var] = str(threads)

    output_dir = ensure_output_dir(cfg.get("output_dir", "outputs"))

    # -------------------------------
    # Load monomer and center it
    # -------------------------------
    input_xyz = cfg.get("input_xyz", "alkyl-pdi.xyz")
    print(f"Loading molecule from {input_xyz} ...")
    monomer_coords, atom_types = read_xyz_file(input_xyz)
    n_atoms = len(atom_types)
    print(f"Loaded {n_atoms} atoms: {', '.join(sorted(set(atom_types)))}")

    monomer_center = np.mean(monomer_coords, axis=0)
    monomer_coords_centered = monomer_coords - monomer_center
    print(f"Centered molecule (original center: {monomer_center})")

    # -------------------------------
    # Configure xTB
    # -------------------------------
    xcfg = XTBConfig(**cfg["xtb"])
    xtb = XTBRunner(xcfg)

    # Precompute monomer energy (centered monomer)
    print("Computing monomer single-point energy with xTB...")
    E_mono_h = xtb.single_point_hartree(atom_types, monomer_coords_centered)
    E_mono_kj = E_mono_h * HARTREE_TO_KJMOL
    print(f"  E_monomer = {E_mono_h: .8f} Eh  ({E_mono_kj: .3f} kJ/mol)")

    # -------------------------------
    # Build the objective function (n-molecule binding energy)
    # -------------------------------
    objective_n = int(cfg.get("objective_n", 2))
    if objective_n < 2:
        raise ValueError(
            "`objective_n` must be >= 2 for a meaningful binding energy.")

    penalty_weight = 2.0
    clash_cutoff = 1.6

    def build_objective_for_pso(pso_cfg: PSOConfig) -> Callable[[np.ndarray], float]:
        # This closure allows (future) dependence on PSO config if desired.
        def objective(vec: np.ndarray) -> float:
            return n_body_binding_energy_kj(
                vec, xtb, monomer_coords_centered, atom_types, E_mono_h,
                n_layers=objective_n,
                penalty_weight=penalty_weight,
                clash_cutoff=clash_cutoff
            )
        return objective

    # Optional parameter bounds
    bounds_cfg = cfg.get("bounds", None)
    bounds: Optional[Dict[str, Tuple[float, float]]] = None
    if isinstance(bounds_cfg, dict):
        bounds = {k: (float(v[0]), float(v[1])) for k, v in bounds_cfg.items()}

    # -------------------------------
    # Hyperparameter tuning (random search)
    # -------------------------------
    # Merge base PSO config + worker count from top-level cfg
    base_pso_cfg = dict(cfg["pso"])
    base_pso_cfg["pso_workers"] = int(cfg.get("pso_workers", 1))

    tuning_cfg = TuningConfig(**cfg.get("tuning", {}))
    if tuning_cfg.enabled:
        print("\nTuning PSO hyperparameters (random search)...")
        tuned_pso_cfg = tune_pso_hyperparameters(
            base_pso_cfg=base_pso_cfg,
            tcfg=tuning_cfg,
            build_objective=build_objective_for_pso,
            bounds=bounds,
            out_dir=output_dir
        )
        print("  Tuning complete. Selected hyperparameters:")
        for k in ["swarm_size", "inertia", "cognitive", "social"]:
            print(f"    {k}: {tuned_pso_cfg[k]}")
    else:
        tuned_pso_cfg = base_pso_cfg

    # -------------------------------
    # Final PSO run with tuned hyperparameters
    # -------------------------------
    pso_cfg = PSOConfig(**tuned_pso_cfg)
    pso = PSO(pso_cfg)

    # Prepare progress loggers (the three required output files)
    energy_log_path = os.path.join(output_dir, "energy_log.csv")
    swarm_trace_path = os.path.join(output_dir, "swarm_trace.csv")
    results_path = os.path.join(output_dir, "optimization_results.txt")

    energy_logger = EnergyLogger(energy_log_path)
    swarm_logger = SwarmTraceLogger(swarm_trace_path)

    # Make the objective callable
    score_fn = build_objective_for_pso(pso.cfg)

    # Define per-iteration hook to log both aggregate and per-particle data.
    def on_iter(it: int, P: np.ndarray, V: np.ndarray, fitness: np.ndarray):
        energy_logger.log(it, fitness)
        swarm_logger.log(it, P, fitness)
        if pso_cfg.verbose_every and (it % pso_cfg.verbose_every == 0 or it == 1):
            best = float(np.min(fitness))
            mean = float(np.mean(fitness))
            std = float(np.std(fitness))
            print(
                f"[Trace] iter={it:4d}  best={best:.6f}  mean={mean:.6f}  std={std:.6f}")

    print(f"\nStarting PSO optimization (xTB objective, n={objective_n}) ...")
    print(f"Configuration: swarm={pso_cfg.swarm_size}, max_iters={pso_cfg.max_iters}, "
          f"workers={pso_cfg.pso_workers}, threads={threads}")

    t0 = time.time()
    best_params, best_val = pso.optimize(
        score_fn, bounds=bounds, on_iter=on_iter)
    t1 = time.time()

    print(f"\n{'='*72}")
    print("OPTIMIZATION RESULTS (xTB, n-body binding energy)")
    print(f"{'='*72}")
    print(f"Best objective (ΔE_bind(n) + penalty): {best_val:.6f} kJ/mol")
    print(f"Best parameters: {summarize_params(best_params)}")
    print(f"Elapsed: {t1 - t0:.1f} s")

    # Persist the results (with the required 'OPTIMIZATION RESULTS' marker at top)
    write_optimization_results(
        results_path, best_params, best_val, pso_cfg, objective_n)

    # -------------------------------
    # Build final stack (using possibly larger `build_n`) and write XYZ
    # -------------------------------
    build_n = int(cfg.get("build_n", 10))
    if build_n < 2:
        print(
            f"Requested build_n={build_n} < 2; forcing to 2 for a meaningful stack.")
        build_n = 2

    best_M = transformation_matrix(best_params)
    # Build layers in *centered* space, then shift back to original placement
    centered_stack_coords, stack_atoms = build_stack_from_transform(
        monomer_coords_centered, atom_types, best_M, n_layers=build_n
    )
    final_stack_coords = centered_stack_coords + monomer_center  # shift back
    total_atoms = final_stack_coords.shape[0]

    # stack statistics (per-layer z-centroids)
    layer_size = monomer_coords.shape[0]
    layer_z = [
        layer_centroid_z(final_stack_coords[i*layer_size:(i+1)*layer_size, :])
        for i in range(build_n)
    ]
    stats = stack_statistics_from_centroids(layer_z)

    print("\nStack Statistics (from per-layer z-centroids):")
    print(f"  Total atoms: {total_atoms}")
    print(f"  Stack height: {stats['stack_height']:.3f} Å")
    print(
        f"  Average layer separation: {stats['avg_sep']:.3f} ± {stats['std_sep']:.3f} Å")

    # Save the complete stack
    xyz_path = os.path.join(
        output_dir, f"molecular_stack_{build_n}molecules.xyz")
    write_xyz_file(
        xyz_path,
        final_stack_coords,
        stack_atoms,
        comment=(f"{build_n}-molecule stack (best objective={best_val:.6f} kJ/mol; "
                 f"avg_sep={stats['avg_sep']:.3f} Å ± {stats['std_sep']:.3f} Å)")
    )
    print(f"\nOutput files generated in {output_dir}:")
    print(f"  • energy_log.csv              (progress log)")
    print(f"  • optimization_results.txt    (best params + tuned PSO hyperparameters)")
    print(f"  • swarm_trace.csv             (per-particle trace at each iteration)")
    print(f"  • {os.path.basename(xyz_path)} (final stacked geometry)")
    print("\nMolecular stack optimization with xTB completed.")


# ============================================================================
# Entry point
# ============================================================================

if __name__ == "__main__":
    import sys
    user_cfg_path = sys.argv[1] if len(sys.argv) > 1 else None
    main(user_cfg_path)
