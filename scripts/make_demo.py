"""Cut a small demo trajectory (protein + ligand, no solvent) out of a full MD run.

Usage: python scripts/make_demo.py TOPOLOGY.prmtop TRAJ.dcd OUT_PREFIX [--stride N]
"""
import argparse
import MDAnalysis as mda
from MDAnalysis import transformations as T

ap = argparse.ArgumentParser()
ap.add_argument("top")
ap.add_argument("traj")
ap.add_argument("out")
ap.add_argument("--stride", type=int, default=2)
ap.add_argument("--keep", default="protein or resname MOL")
a = ap.parse_args()

u = mda.Universe(a.top, a.traj)
keep = u.select_atoms(a.keep)
prot = u.select_atoms("protein")
u.trajectory.add_transformations(
    T.unwrap(u.atoms), T.center_in_box(prot, wrap=False), T.wrap(keep, compound="fragments")
)
keep.write(f"{a.out}.pdb", frames=u.trajectory[:1])
with mda.Writer(f"{a.out}.xtc", keep.n_atoms) as w:
    for _ in u.trajectory[:: a.stride]:
        w.write(keep)
print(f"{keep.n_atoms} atoms, {len(u.trajectory[:: a.stride])} frames -> {a.out}.pdb/.xtc")
