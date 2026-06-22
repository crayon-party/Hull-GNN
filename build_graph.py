import networkx as nx
import matplotlib.pyplot as plt
import json

from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCC.Core.TopTools import (
    TopTools_IndexedMapOfShape,
    TopTools_IndexedDataMapOfShapeListOfShape,
)
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop

STEP_FILE = "GPPH-A1-B3-Hull.step"   # <-- change this

reader = STEPControl_Reader()
assert reader.ReadFile(STEP_FILE) == IFSelect_RetDone
reader.TransferRoots()
shape = reader.OneShape()

# --- 1. Collect the faces (these are our nodes) ---
faces = []
exp = TopExp_Explorer(shape, TopAbs_FACE)
while exp.More():
    faces.append(exp.Current())
    exp.Next()
print(f"Collected {len(faces)} faces (nodes)")

# --- 2. Per-face area (first node feature) ---
def face_area(face):
    props = GProp_GProps()
    brepgprop.SurfaceProperties(face, props)
    return props.Mass()

for i, f in enumerate(faces):
    print(f"  node {i}: area = {face_area(f)/1e6:.3f} m^2")

# --- 3. Build adjacency: two faces are connected if they share an edge ---
# Map each edge -> the list of faces using it
# --- 3. Build adjacency: two faces share an edge ---
from OCC.Core.TopTools import (
    TopTools_IndexedDataMapOfShapeListOfShape,
    TopTools_IndexedMapOfShape,
    TopTools_ListIteratorOfListOfShape,
)
from OCC.Core.TopExp import topexp

edge_face_map = TopTools_IndexedDataMapOfShapeListOfShape()
topexp.MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edge_face_map)

face_map = TopTools_IndexedMapOfShape()
for f in faces:
    face_map.Add(f)

edges_set = set()
shared_edge_lengths = {}

n_edges = edge_face_map.Size()          # <-- was .Extent()
for i in range(1, n_edges + 1):
    edge = edge_face_map.FindKey(i)
    adj_faces = edge_face_map.FindFromIndex(i)
    idxs = []
    lit = TopTools_ListIteratorOfListOfShape(adj_faces)
    while lit.More():
        fidx = face_map.FindIndex(lit.Value())
        if fidx > 0:
            idxs.append(fidx - 1)        # 0-based node index
        lit.Next()
    if len(idxs) == 2:
        a, b = sorted(idxs)
        props = GProp_GProps()
        brepgprop.LinearProperties(edge, props)
        key = (a, b)
        edges_set.add(key)
        shared_edge_lengths[key] = shared_edge_lengths.get(key, 0.0) + props.Mass()

print(f"\nFound {len(edges_set)} graph edges (seams between faces):")
for (a, b) in sorted(edges_set):
    print(f"  {a} <-> {b}   shared boundary length = {shared_edge_lengths[(a,b)]/1e3:.2f} m")


# --- Build a networkx graph from what we extracted ---
G = nx.Graph()

# Nodes: one per face, with area as a feature
for i, f in enumerate(faces):
    G.add_node(i, area=face_area(f) / 1e6)  # area in m^2

# Edges: one per shared boundary, with weld_length as a feature
for (a, b), length in shared_edge_lengths.items():
    G.add_edge(a, b, weld_length=length / 1e3)  # length in m

print(f"\nGraph built: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
print("Connected:", nx.is_connected(G))
print("Node degrees:", dict(G.degree()))

# --- Save it so you never have to re-run the geometry ---
# (1) as JSON — human-readable, portable
data = nx.node_link_data(G)
with open("hull_graph.json", "w") as fp:
    json.dump(data, fp, indent=2)
print("Saved hull_graph.json")

# (2) as GraphML — opens in tools like yEd/Gephi if you want
nx.write_graphml(G, "hull_graph.graphml")

# --- Draw it ---
#pos = nx.spring_layout(G, seed=42, k=0.8)  # auto-layout

# then pass node_color=colors to draw_networkx_nodes
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop

def face_centroid(face):
    props = GProp_GProps()
    brepgprop.SurfaceProperties(face, props)
    c = props.CentreOfMass()
    return (c.X(), c.Y(), c.Z())

# Real 3D centroids, then project to 2D for plotting.
# A hull is long (X) and tall (Z); side view = (X, Z) usually reads best.
centroids = {i: face_centroid(f) for i, f in enumerate(faces)}
pos = {i: (c[0], c[2]) for i, c in centroids.items()}   # (X, Z) side view
areas = [G.nodes[n]["area"] for n in G.nodes]
sizes = [80 + a * 60 for a in areas]  # node size ~ area
weights = [G.edges[e]["weld_length"] for e in G.edges]

plt.figure(figsize=(10, 8))
nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color="#5DCAA5")
nx.draw_networkx_edges(G, pos, width=[0.5 + w * 0.2 for w in weights], alpha=0.5)
nx.draw_networkx_labels(G, pos, font_size=10)
plt.title("Hull face-adjacency graph")
plt.axis("off")
plt.tight_layout()
plt.savefig("hull_graph.png", dpi=150)
print("Saved hull_graph.png")
plt.show()


pos_spring = nx.spring_layout(G, seed=42, k=1.5, iterations=200)

plt.figure(figsize=(10, 8))
nx.draw_networkx_nodes(G, pos_spring, node_size=sizes, node_color="#5DCAA5")
nx.draw_networkx_edges(G, pos_spring, width=[0.5 + w * 0.2 for w in weights], alpha=0.5)
nx.draw_networkx_labels(G, pos_spring, font_size=9)
plt.title("Hull face-adjacency graph (spring layout)")
plt.axis("off")
plt.tight_layout()
plt.savefig("hull_graph_spring.png", dpi=150)
plt.show()



# =====================================================================
#  CURVATURE EXTRACTION
#  --------------------------------------------------------------------
#  For each face we sample a grid of points across its (u,v) surface
#  and measure the two PRINCIPAL CURVATURES at each point.
#    - Gaussian curvature K = k1 * k2   (product)
#    - Mean curvature     H = (k1+k2)/2 (average)
#
#  Why these two:
#    K tells us DEVELOPABILITY. K≈0 everywhere => developable (flat or
#    single-curvature, cheap to form). K far from 0 => double curvature
#    (needs line heating, expensive). K changing SIGN across the face
#    => reverse double curvature (saddle<->dome, the worst case).
#    H tells us the overall BENDING MAGNITUDE.
#
#  We will use the per-face SUMMARY of these (mean, max, sign-change
#  fraction) two ways:
#    (1) ONE-TIME SANITY CHECK  -> confirm extraction is correct
#    (2) KEPT NODE FEATURE       -> describes plate shape for the model
# =====================================================================

from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepLProp import BRepLProp_SLProps
from OCC.Core.TopAbs import TopAbs_REVERSED
import numpy as np

def face_curvature_summary(face, n_samples=8):
    """
    Sample an n_samples x n_samples grid over the face's (u,v) domain,
    compute principal curvatures at each valid point, and return a
    summary dict of Gaussian (K) and mean (H) curvature statistics.
    """
    # BRepAdaptor_Surface wraps the face so we can query its (u,v) range
    # and evaluate geometry/curvature at any (u,v) parameter.
    surf = BRepAdaptor_Surface(face)

    # Parameter bounds of this face in u and v directions.
    u_min, u_max = surf.FirstUParameter(), surf.LastUParameter()
    v_min, v_max = surf.FirstVParameter(), surf.LastVParameter()

    # SLProps computes Surface Local PROPerties (curvature, normals...).
    # Args: (surface, n_derivatives=2, resolution/tolerance).
    # We need 2 derivatives because curvature is a 2nd-order quantity.
    props = BRepLProp_SLProps(surf, 2, 1e-6)

    K_vals = []   # Gaussian curvature samples
    H_vals = []   # Mean curvature samples

    # Walk a regular grid across the (u,v) domain.
    for i in range(n_samples):
        for j in range(n_samples):
            # Map grid index -> actual (u,v) parameter on the face.
            # The +1 spacing keeps us just inside the edges, avoiding
            # the boundary where curvature can be ill-defined.
            u = u_min + (u_max - u_min) * (i + 0.5) / n_samples
            v = v_min + (v_max - v_min) * (j + 0.5) / n_samples

            props.SetParameters(u, v)

            # Curvature is only defined where the surface is smooth and
            # the local frame is well-defined. Skip points where it isn't.
            if not props.IsCurvatureDefined():
                continue

            # OpenCASCADE gives the two principal curvatures directly.
            k1 = props.MaxCurvature()
            k2 = props.MinCurvature()

            # NOTE on sign: a face can be "reversed" relative to its
            # surface (the TopAbs_REVERSED orientation). That flips the
            # normal direction and hence the SIGN of curvature. We flip
            # it back so all faces share a consistent sign convention --
            # important later when sign of K marks reverse curvature.
            if face.Orientation() == TopAbs_REVERSED:
                k1, k2 = -k1, -k2

            K_vals.append(k1 * k2)        # Gaussian
            H_vals.append((k1 + k2) / 2)  # Mean

    # If no valid samples (rare: tiny sliver faces), return zeros so the
    # pipeline doesn't crash -- we'll see these as 0 and can inspect them.
    if len(K_vals) == 0:
        return dict(K_mean=0.0, K_absmean=0.0, K_max=0.0,
                    H_absmean=0.0, sign_change_frac=0.0, n_valid=0)

    K_vals = np.array(K_vals)
    H_vals = np.array(H_vals)

    # sign_change_frac: fraction of samples whose Gaussian curvature sign
    # differs from the face's dominant sign. High value => the face
    # contains BOTH saddle and dome regions => reverse double curvature.
    dominant_sign = np.sign(np.median(K_vals)) or 1
    sign_change_frac = np.mean(np.sign(K_vals) != dominant_sign)

    return dict(
        K_mean=float(K_vals.mean()),          # signed average Gaussian
        K_absmean=float(np.abs(K_vals).mean()),  # how non-developable overall
        K_max=float(np.abs(K_vals).max()),    # worst local double-curvature
        H_absmean=float(np.abs(H_vals).mean()),  # overall bending magnitude
        sign_change_frac=float(sign_change_frac),  # reverse-curvature flag
        n_valid=int(len(K_vals)),             # how many samples were usable
    )

# Compute the summary for every face.
curv = {i: face_curvature_summary(f) for i, f in enumerate(faces)}

# ---------------------------------------------------------------------
#  (1) ONE-TIME SANITY CHECK  -- throwaway, just to trust the extraction
#  Expectation from naval architecture:
#    big side panels  -> LOW  K_absmean (near-developable midship shell)
#    bow/stern faces  -> HIGH K_absmean and/or sign changes (compound)
# ---------------------------------------------------------------------
print("\n=== Curvature sanity check (units: 1/mm and 1/mm^2) ===")
print(f"{'node':>4} {'area(m2)':>9} {'K_absmean':>12} {'K_max':>12} "
      f"{'H_absmean':>10} {'signflip':>9} {'nval':>5}")
for i in range(len(faces)):
    a = face_area(faces[i]) / 1e6
    c = curv[i]
    print(f"{i:>4} {a:>9.3f} {c['K_absmean']:>12.3e} {c['K_max']:>12.3e} "
          f"{c['H_absmean']:>10.3e} {c['sign_change_frac']:>9.2f} {c['n_valid']:>5}")

# ---------------------------------------------------------------------
#  (2) KEPT NODE FEATURE -- attach curvature to each node in G.
#  (Make sure this runs AFTER G is created and nodes are added.)
# ---------------------------------------------------------------------
for i in range(len(faces)):
    c = curv[i]
    G.nodes[i]["K_absmean"]        = c["K_absmean"]
    G.nodes[i]["K_max"]            = c["K_max"]
    G.nodes[i]["H_absmean"]        = c["H_absmean"]
    G.nodes[i]["sign_change_frac"] = c["sign_change_frac"]

print("\nCurvature features attached to graph nodes.")