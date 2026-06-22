# surface_test.py
# 1) Reports the surface TYPE of each face (Plane / BSpline / etc.)
# 2) Probes faces whose curvature evaluates to zero, comparing them
#    against a known-good curved face, to see why they fail.

from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_REVERSED
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepLProp import BRepLProp_SLProps
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.GeomAbs import (
    GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Sphere,
    GeomAbs_Torus, GeomAbs_BezierSurface, GeomAbs_BSplineSurface,
    GeomAbs_SurfaceOfRevolution, GeomAbs_SurfaceOfExtrusion,
    GeomAbs_OffsetSurface, GeomAbs_OtherSurface,
)
import numpy as np

STEP_FILE = "GPPH-A1-B3-Hull.step"   # <-- change this

# --- Load ---
reader = STEPControl_Reader()
assert reader.ReadFile(STEP_FILE) == IFSelect_RetDone, "Failed to read STEP file"
reader.TransferRoots()
shape = reader.OneShape()

faces = []
exp = TopExp_Explorer(shape, TopAbs_FACE)
while exp.More():
    faces.append(exp.Current()); exp.Next()

def face_area(face):
    p = GProp_GProps(); brepgprop.SurfaceProperties(face, p); return p.Mass()

type_names = {
    GeomAbs_Plane: "Plane", GeomAbs_Cylinder: "Cylinder",
    GeomAbs_Cone: "Cone", GeomAbs_Sphere: "Sphere",
    GeomAbs_Torus: "Torus", GeomAbs_BezierSurface: "Bezier",
    GeomAbs_BSplineSurface: "BSpline",
    GeomAbs_SurfaceOfRevolution: "Revolution",
    GeomAbs_SurfaceOfExtrusion: "Extrusion",
    GeomAbs_OffsetSurface: "Offset", GeomAbs_OtherSurface: "Other",
}

# ============================================================
#  PART 1 — surface type per face + mean curvature (to flag zeros)
# ============================================================
def mean_curv(face, n=5):
    """Return mean |k1|+|k2| magnitude over a grid; ~0 means it failed."""
    surf = BRepAdaptor_Surface(face)
    u0, u1 = surf.FirstUParameter(), surf.LastUParameter()
    v0, v1 = surf.FirstVParameter(), surf.LastVParameter()
    props = BRepLProp_SLProps(surf, 2, 1e-9)
    mag = []
    for i in range(n):
        for j in range(n):
            u = u0 + (u1-u0)*(i+0.5)/n
            v = v0 + (v1-v0)*(j+0.5)/n
            props.SetParameters(u, v)
            if props.IsCurvatureDefined():
                mag.append(abs(props.MaxCurvature()) + abs(props.MinCurvature()))
    return np.mean(mag) if mag else 0.0

print(f"\n{'node':>4}  {'type':<10} {'area(m2)':>9} {'curv_mag':>11}")
print("-" * 40)
zero_faces, good_faces = [], []
for i, f in enumerate(faces):
    stype = type_names.get(BRepAdaptor_Surface(f).GetType(), "?")
    cm = mean_curv(f)
    print(f"{i:>4}  {stype:<10} {face_area(f)/1e6:>9.3f} {cm:>11.3e}")
    # flag big-ish faces whose curvature came back essentially zero
    if cm < 1e-12 and face_area(f)/1e6 > 1.0:
        zero_faces.append(i)
    elif cm > 1e-6:
        good_faces.append(i)

# ============================================================
#  PART 2 — diagnose the zero faces vs a known-good face
# ============================================================
if zero_faces:
    ref = good_faces[0] if good_faces else None
    probe_list = ([ref] if ref is not None else []) + zero_faces
    print(f"\n=== Zero-curvature diagnostic ===")
    print(f"(comparing zero faces {zero_faces} against good face {ref})")
    for idx in probe_list:
        f = faces[idx]
        surf = BRepAdaptor_Surface(f)
        u0, u1 = surf.FirstUParameter(), surf.LastUParameter()
        v0, v1 = surf.FirstVParameter(), surf.LastVParameter()
        props = BRepLProp_SLProps(surf, 2, 1e-9)
        props.SetParameters((u0+u1)/2, (v0+v1)/2)
        tag = "GOOD" if idx == ref else "ZERO"
        print(f"\n--- face {idx} [{tag}] ---")
        print(f"  u: {u0:.4f}..{u1:.4f}   v: {v0:.4f}..{v1:.4f}")
        print(f"  curvature defined at center: {props.IsCurvatureDefined()}")
        if props.IsCurvatureDefined():
            print(f"  k_max: {props.MaxCurvature():.3e}  k_min: {props.MinCurvature():.3e}")
        print(f"  tangentU defined: {props.IsTangentUDefined()}")
else:
    print("\nNo zero-curvature faces flagged — curvature looks usable on all sizable faces.")
