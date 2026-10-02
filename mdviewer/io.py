"""Loading and rendering helpers that sit between MDAnalysis and the UI."""
from __future__ import annotations

import tempfile
from pathlib import Path

import MDAnalysis as mda

TOPOLOGY_TYPES = ["pdb", "gro", "prmtop", "parm7", "psf"]
TRAJECTORY_TYPES = ["xtc", "dcd", "trr", "nc"]
MAX_ANIMATION_FRAMES = 100


def save_upload(uploaded, dirpath: Path) -> Path:
    path = dirpath / Path(uploaded.name).name
    path.write_bytes(uploaded.getbuffer())
    return path


def load_universe(topology: str | Path, trajectory: str | Path) -> mda.Universe:
    return mda.Universe(str(topology), str(trajectory))


def new_workdir() -> Path:
    return Path(tempfile.mkdtemp(prefix="mdviewer-"))


def frame_pdb(u: mda.Universe, frame: int, select: str = "all") -> str:
    """PDB text of one frame for the selected atoms."""
    ag = u.select_atoms(select)
    u.trajectory[frame]
    return _write_pdb(ag, [frame], u)


def multi_model_pdb(u: mda.Universe, select: str = "all", stride: int = 1) -> str:
    """Multi-model PDB for browser-side playback, capped at MAX_ANIMATION_FRAMES."""
    ag = u.select_atoms(select)
    n = len(u.trajectory)
    step = max(stride, -(-n // MAX_ANIMATION_FRAMES))
    return _write_pdb(ag, range(0, n, step), u)


def _write_pdb(ag, frames, u) -> str:
    # MDAnalysis' PDB writer closes the file object it is given, so write to disk and read back.
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "frames.pdb"
        with mda.Writer(str(path), ag.n_atoms, multiframe=True) as w:
            for i in frames:
                u.trajectory[i]
                w.write(ag)
        return path.read_text()
