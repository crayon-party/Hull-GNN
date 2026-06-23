from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop

STEP_FILE = "GPPH-A1-B3-Hull.step"   # <-- change this

reader = STEPControl_Reader()
assert reader.ReadFile(STEP_FILE) == IFSelect_RetDone
reader.TransferRoots()
shape = reader.OneShape()

faces = []
exp = TopExp_Explorer(shape, TopAbs_FACE)
while exp.More():
    faces.append(exp.Current()); exp.Next()

def centroid(face):
    p = GProp_GProps()
    brepgprop.SurfaceProperties(face, p)
    c = p.CentreOfMass()
    return (c.X(), c.Y(), c.Z())

def face_area(face):
    p = GProp_GProps()
    brepgprop.SurfaceProperties(face, p)
    return p.Mass()

print(f"{'node':>4} {'area(m2)':>9} {'cx':>10} {'cy':>10} {'cz':>10}")
print("-" * 54)
for i, f in enumerate(faces):
    cx, cy, cz = centroid(f)
    print(f"{i:>4} {face_area(f)/1e6:>9.3f} {cx:>10.1f} {cy:>10.1f} {cz:>10.1f}")