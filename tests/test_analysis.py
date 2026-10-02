import numpy as np
import MDAnalysis as mda
import pytest
from MDAnalysis.coordinates.memory import MemoryReader

from mdviewer import analysis, io


def make_universe(coords):
    n_frames, n_atoms, _ = coords.shape
    u = mda.Universe.empty(n_atoms, n_residues=n_atoms, atom_resindex=np.arange(n_atoms),
                           trajectory=True)
    u.add_TopologyAttr("name", ["CA"] * n_atoms)
    u.add_TopologyAttr("resname", ["ALA"] * n_atoms)
    u.add_TopologyAttr("resid", np.arange(1, n_atoms + 1))
    u.add_TopologyAttr("segid", ["A"])
    u.add_TopologyAttr("masses", np.ones(n_atoms))
    u.load_new(coords.astype(np.float32), format=MemoryReader)
    return u


@pytest.fixture
def rigid():
    base = np.random.default_rng(0).normal(size=(10, 3)) * 5
    return make_universe(np.stack([base] * 6))


def test_frame_indices():
    assert list(analysis.frame_indices(10, 3)) == [0, 3, 6, 9]
    with pytest.raises(ValueError):
        analysis.frame_indices(10, 0)


def test_rigid_trajectory_has_zero_rmsd_and_rmsf(rigid):
    assert np.allclose(analysis.rmsd_series(rigid, "name CA"), 0, atol=1e-4)
    _, f = analysis.rmsf_per_residue(rigid, "name CA")
    assert np.allclose(f, 0, atol=1e-4)


def test_rmsd_detects_displacement():
    base = np.random.default_rng(1).normal(size=(10, 3)) * 5
    moved = base.copy()
    moved[0] += [3, 0, 0]
    r = analysis.rmsd_series(make_universe(np.stack([base, moved])), "name CA")
    assert r[0] == pytest.approx(0, abs=1e-4) and r[1] > 0.5


def test_rg_known_value():
    pts = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], float)
    u = make_universe(np.stack([pts, pts * 2]))
    assert analysis.rg_series(u, "all") == pytest.approx([1.0, 2.0], abs=1e-4)


def test_stride_reduces_frames(rigid):
    assert len(analysis.rg_series(rigid, "all", stride=2)) == 3


def test_empty_selection_raises(rigid):
    for fn in (analysis.rmsd_series, analysis.rg_series, analysis.rmsf_per_residue):
        with pytest.raises(ValueError):
            fn(rigid, "resname XXX")


def test_multi_model_pdb_caps_frames():
    pts = np.random.default_rng(2).normal(size=(250, 4, 3))
    pdb = io.multi_model_pdb(make_universe(pts))
    assert pdb.count("\nMODEL") + pdb.startswith("MODEL") <= io.MAX_ANIMATION_FRAMES


def test_demo_data_loads():
    from pathlib import Path
    d = Path(__file__).parent.parent / "demo"
    u = io.load_universe(d / "egfr_T790M_osimertinib.pdb", d / "egfr_T790M_osimertinib.xtc")
    assert len(u.trajectory) == 250 and len(u.select_atoms("resname MOL")) > 0
