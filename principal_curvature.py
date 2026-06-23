# principal_curvature.py
# TRUE principal curvatures via the first/second fundamental forms,
# evaluated from 3D surface derivatives (faithful to Nam's method),
# computed directly so it doesn't depend on (u,v) parameter scaling.

from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_REVERSED
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.gp import gp_Pnt, gp_Vec
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

def face_area(face):
    p = GProp_GProps(); brepgprop.SurfaceProperties(face, p); return p.Mass()

def principal_curvatures_at(surface, u, v):
    """
    Compute (k1, k2) at parameter (u,v) from the fundamental forms.
    Uses D2: position + 1st and 2nd derivatives. Returns None if degenerate.
    """
    p   = gp_Pnt()
    du  = gp_Vec(); dv  = gp_Vec()      # first derivatives r_u, r_v
    duu = gp_Vec(); dvv = gp_Vec(); duv = gp_Vec()  # second derivatives
    surface.D2(u, v, p, du, dv, duu, duv, dvv)

    ru  = np.array([du.X(),  du.Y(),  du.Z()])
    rv  = np.array([dv.X(),  dv.Y(),  dv.Z()])
    ruu = np.array([duu.X(), duu.Y(), duu.Z()])
    ruv = np.array([duv.X(), duv.Y(), duv.Z()])
    rvv = np.array([dvv.X(), dvv.Y(), dvv.Z()])

    # Unit normal from the cross product of the tangents.
    n = np.cross(ru, rv)
    nlen = np.linalg.norm(n)
    if nlen < 1e-12:
        return None          # degenerate point (tangents parallel)
    n /= nlen

    # First fundamental form
    E = ru @ ru; F = ru @ rv; G = rv @ rv
    # Second fundamental form
    L = ruu @ n; M = ruv @ n; N = rvv @ n

    denom = E*G - F*F
    if abs(denom) < 1e-20:
        return None
    K = (L*N - M*M) / denom                 # Gaussian
    H = (E*N - 2*F*M + G*L) / (2*denom)      # Mean
    disc = max(H*H - K, 0.0)                 # clamp tiny negatives
    root = np.sqrt(disc)
    k1, k2 = H + root, H - root
    return k1, k2

def face_curvature(face, n=8):
    surf = BRep_Tool.Surface(face)          # underlying Geom_Surface
    # parameter bounds
    umin, umax, vmin, vmax = BRep_Tool.Surface(face).Bounds() \
        if hasattr(BRep_Tool.Surface(face), "Bounds") else (0,1,0,1)
    # safer bounds via the face:
    from OCC.Core.BRepTools import breptools
    umin, umax, vmin, vmax = breptools.UVBounds(face)

    reversed_ = (face.Orientation() == TopAbs_REVERSED)
    k1s, k2s = [], []
    for i in range(n):
        for j in range(n):
            u = umin + (umax-umin)*(i+0.5)/n
            v = vmin + (vmax-vmin)*(j+0.5)/n
            res = principal_curvatures_at(surf, u, v)
            if res is None:
                continue
            k1, k2 = res
            if reversed_:
                k1, k2 = -k1, -k2
            k1s.append(k1); k2s.append(k2)
    if not k1s:
        return None
    k1s = np.array(k1s); k2s = np.array(k2s)
    return dict(
        n_valid=len(k1s),
        k1_mean=k1s.mean(), k2_mean=k2s.mean(),
        K_mean=(k1s*k2s).mean(),              # mean Gaussian
        H_mean=((k1s+k2s)/2).mean(),          # mean mean-curvature
        kmax_abs=np.abs(np.concatenate([k1s,k2s])).max(),
    )

print(f"{'node':>4} {'area':>8} {'nval':>5} "
      f"{'k1_mean':>11} {'k2_mean':>11} {'R1(m)':>9} {'R2(m)':>9}")
print("-"*64)
for i, f in enumerate(faces):
    c = face_curvature(f)
    a = face_area(f)/1e6
    if c is None:
        print(f"{i:>4} {a:>8.3f}   -- no valid samples --")
        continue
    R1 = 1.0/abs(c['k1_mean'])/1000 if abs(c['k1_mean'])>1e-30 else float('inf')
    R2 = 1.0/abs(c['k2_mean'])/1000 if abs(c['k2_mean'])>1e-30 else float('inf')
    print(f"{i:>4} {a:>8.3f} {c['n_valid']:>5} "
          f"{c['k1_mean']:>11.3e} {c['k2_mean']:>11.3e} {R1:>9.1f} {R2:>9.1f}")