"""Pure analysis helpers: no Streamlit, so they can be unit tested."""
from __future__ import annotations

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms


def frame_indices(n_frames: int, stride: int = 1) -> np.ndarray:
    if stride < 1:
        raise ValueError("stride must be >= 1")
    return np.arange(0, n_frames, stride)


def rmsd_series(u: mda.Universe, select: str = "protein and name CA",
                align_on: str | None = None, stride: int = 1) -> np.ndarray:
    """RMSD in Å of `select` against frame 0, after superposing on `align_on` (default: select)."""
    ag = u.select_atoms(select)
    ref = u.select_atoms(align_on or select)
    if len(ag) == 0 or len(ref) == 0:
        raise ValueError(f"selection matched no atoms: {select!r}")
    runner = rms.RMSD(ag, ref=ag, select=align_on or select, groupselections=None)
    runner.run(step=stride)
    return runner.results.rmsd[:, 2]


def rg_series(u: mda.Universe, select: str = "protein", stride: int = 1) -> np.ndarray:
    """Radius of gyration in Å over the trajectory."""
    ag = u.select_atoms(select)
    if len(ag) == 0:
        raise ValueError(f"selection matched no atoms: {select!r}")
    return np.array([ag.radius_of_gyration() for _ in u.trajectory[::stride]])


def rmsf_per_residue(u: mda.Universe, select: str = "protein and name CA",
                     stride: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Per-residue RMSF in Å; trajectory is first fitted on the same selection."""
    ag = u.select_atoms(select)
    if len(ag) == 0:
        raise ValueError(f"selection matched no atoms: {select!r}")
    aligner = align.AlignTraj(u, u, select=select, in_memory=True, step=stride).run()
    del aligner
    r = rms.RMSF(ag).run()
    return ag.resids, r.results.rmsf
