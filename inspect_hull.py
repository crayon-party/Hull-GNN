from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import (
    TopAbs_SOLID, TopAbs_SHELL, TopAbs_FACE, TopAbs_EDGE
)
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop

STEP_FILE = "GPPH-A1-B3-Hull.step"

reader = STEPControl_Reader()
status = reader.ReadFile(STEP_FILE)
assert status == IFSelect_RetDone, "Failed to read STEP file"

reader.TransferRoots()
shape = reader.OneShape()

def count(shape, enum):
    exp = TopExp_Explorer(shape, enum)
    n = 0
    while exp.More():
        n += 1
        exp.Next()
    return n

print("=== Top-level counts ===")
print("Solids:", count(shape, TopAbs_SOLID))
print("Shells:", count(shape, TopAbs_SHELL))
print("Faces :", count(shape, TopAbs_FACE))
print("Edges :", count(shape, TopAbs_EDGE))

# How many faces does each shell contain? (tells us if shells = plates)
print("\n=== Faces per shell (first 20 shells) ===")
exp = TopExp_Explorer(shape, TopAbs_SHELL)
i = 0
while exp.More() and i < 20:
    shell = exp.Current()
    print(f"  shell {i}: {count(shell, TopAbs_FACE)} faces")
    i += 1
    exp.Next()

# Surface area of the whole thing (sanity check on units)
props = GProp_GProps()
brepgprop.SurfaceProperties(shape, props)
print("\n=== Sanity check ===")
print(f"Total surface area: {props.Mass():.2f} (units^2 — tells us mm vs m)")
