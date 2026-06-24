# flatness_test.py
# Independent geometric check: is a face actually flat in 3D?
# Samples real 3D points and measures how far they deviate from the
# best-fit plane. No curvature math -- just "are these points planar?"

from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepTools import breptools
from OCC.Core.gp import gp_Pnt
import numpy as np

STEP_FILE = "GPPH-A1-B3-Hull.step"   # <-- change this

reader = STEPControl_Reader()
assert reader.ReadFile(STEP_FILE) == IFSelect_RetDone
reader.TransferRoots()
shape = reader.OneShape()

faces = []
exp = TopExp_Explorer(shape, TopAbs_FACE)
while exp.More():
    faces.append(exp.Current()); exp.Next()

def sample_points(face, n=10):
    """Return an array of real 3D points sampled across the face."""
    surf = BRep_Tool.Surface(face)
    umin, umax, vmin, vmax = breptools.UVBounds(face)
    pts = []
    for i in range(n):
        for j in range(n):
            u = umin + (umax-umin)*(i+0.5)/n
            v = vmin + (vmax-vmin)*(j+0.5)/n
            p = gp_Pnt()
            surf.D0(u, v, p)           # D0 = just the point, no derivatives
            pts.append([p.X(), p.Y(), p.Z()])
    return np.array(pts)

def planarity(pts):
    """
    Fit the best plane through the points, return max deviation from it.
    ~0 deviation => genuinely flat. Large deviation => curved.
    """
    centroid = pts.mean(axis=0)
    centered = pts - centroid
    # The plane normal is the smallest-singular-value direction (SVD).
    _, s, vh = np.linalg.svd(centered)
    normal = vh[2]                     # least-variance direction = plane normal
    # Distance of each point from the plane:
    dists = np.abs(centered @ normal)
    return dists.max(), dists.mean()

print(f"{'node':>4} {'span(m)':>9} {'maxdev(mm)':>11} {'meandev(mm)':>12}  verdict")
print("-"*56)
for i, f in enumerate(faces):
    pts = sample_points(f)
    span = np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)) / 1000  # m
    maxdev, meandev = planarity(pts)
    # If points deviate from a plane by more than ~1mm, it's curved.
    verdict = "FLAT" if maxdev < 1.0 else "CURVED"
    print(f"{i:>4} {span:>9.2f} {maxdev:>11.2f} {meandev:>12.2f}  {verdict}")