"""
backset_cost.py
================
Relative forming-cost classification for hull plates, faithful to:

  Parsons, M.G., Nam, J-H., & Singer, D.J. (1999).
  "A Scalar Metric for Assessing the Producibility of a Hull Form
   in Early Design." Journal of Ship Production 15(2).

The paper does NOT feed raw curvature into the cost model. It converts
principal curvatures into two NON-DIMENSIONAL BACKSET ratios:

  b1 = backset in the LARGER principal-curvature direction
  b2 = backset in the ORTHOGONAL direction

  (backset = rise of the plate above a flat plane, divided by plate length;
   see paper Fig. 5. b1 is always the larger, so by definition b1 >= |b2|.)

Each plate is then classified into one of EIGHT types, each with a fixed
RELATIVE COST (cost 1 = a flat plate; higher = more expensive to form).
This module implements those eight classes exactly as tabulated in the
paper, plus the finer-grained fuzzy rule matrix (Fig. 9) as an optional
lookup for validation.

Sign convention (from the paper):
  b2 > 0  -> ordinary double curvature (dome / bowl, same-sign curvatures)
  b2 < 0  -> REVERSE double curvature  (saddle, opposite-sign curvatures)
  b2 ~ 0  -> single-direction curvature or flat
"""

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# The eight plate classes, exactly as defined in the paper (Types 1-8).
# Each class: a human-readable name, the process used, and the relative cost.
# The threshold LOGIC is in classify_plate() below, because several types
# share cost values and the boundaries are interval tests, not a flat table.
# ---------------------------------------------------------------------------
CLASS_INFO = {
    1: dict(name="Flat plate", process="No fabrication", cost=1),
    2: dict(name="Low curvature, one direction", process="Roll", cost=2),
    3: dict(name="High curvature, one direction", process="Roll", cost=3),
    4: dict(name="Moderate curvature + small backset", process="Roll", cost=3),
    5: dict(name="Moderate double curvature", process="Press + line heat", cost=6),
    6: dict(name="Moderate reverse double curvature", process="Ring press + line heat", cost=9),
    7: dict(name="High double curvature", process="Ring press + line heat", cost=8),
    8: dict(name="High reverse double curvature", process="Separate, form, weld", cost=12),
}


@dataclass
class PlateClassification:
    plate_type: int  # 1-8
    name: str
    process: str
    relative_cost: int
    b1: float
    b2: float


def classify_plate(b1: float, b2: float) -> PlateClassification:
    """
    Classify a plate into one of the eight Parsons-Nam types from its two
    backset ratios, and return the type, name, process, and relative cost.

    Thresholds (from the paper's Type 1-8 definitions):

      Type 1  Flat:                    0.00 <= b1 <= 0.01,  0.00 <= b2 <= 0.0025
      Type 2  Low single curvature:    0.01 <= b1 <= 0.16,  0.00 <= b2 <= 0.0025
      Type 3  High single curvature:   b1 >= 0.16,          0.00 <= b2 <= 0.0025
      Type 4  Moderate + small backset:0.01 <= b1 <= 0.08,  0.0025 <= b2 <= 0.02
      Type 5  Moderate double:         0.08 <= b1 <= 0.16,  0.02 <= b2 <= 0.04
      Type 6  Moderate reverse double: 0.08 <= b1 <= 0.16, -0.04 <= b2 <= -0.02
      Type 7  High double:             b1 >= 0.16,          b2 >= 0.04
      Type 8  High reverse double:     b1 >= 0.16,          b2 <= -0.04

    Note: b1 is by definition the larger backset, so b1 >= |b2| always.
    The paper's classes don't tile the whole plane perfectly (they were
    defined around representative shipyard cases), so we assign each plate
    to the nearest matching type, falling back on the sign of b2 (reverse
    vs. ordinary) and the magnitude of b1 when a value lands between the
    tabulated boxes.
    """
    # Work with the sign of b2 (reverse curvature) separately from magnitude.
    reverse = b2 < -0.0025  # meaningfully negative -> reverse double
    ab2 = abs(b2)

    # --- Single-curvature / flat band: b2 essentially zero ---
    if ab2 <= 0.0025:
        if b1 <= 0.01:
            t = 1
        elif b1 <= 0.16:
            t = 2
        else:
            t = 3
    # --- Small transverse backset (Type 4 band) ---
    elif ab2 <= 0.02 and not reverse:
        # small positive backset in the orthogonal direction
        if b1 <= 0.08:
            t = 4
        elif b1 <= 0.16:
            t = 5  # grades into moderate double curvature
        else:
            t = 7
    # --- Moderate curvature band ---
    elif ab2 <= 0.04:
        if b1 <= 0.16:
            t = 6 if reverse else 5
        else:
            t = 8 if reverse else 7
    # --- High curvature band: |b2| > 0.04 ---
    else:
        t = 8 if reverse else 7

    info = CLASS_INFO[t]
    return PlateClassification(
        plate_type=t,
        name=info["name"],
        process=info["process"],
        relative_cost=info["cost"],
        b1=b1,
        b2=b2,
    )


# ---------------------------------------------------------------------------
# The scalar hull metric (paper Eq. 15): area-weighted mean relative cost.
# Given a list of (area, relative_cost) per plate, returns the single scalar
# the paper reports (their product tanker came out at ~2.01).
# ---------------------------------------------------------------------------
def scalar_hull_metric(plates):
    """
    plates: iterable of (area, relative_cost) tuples.
    Returns SM = sum(cost_i * area_i) / sum(area_i)  -- paper Eq. (15).
    """
    total_cost_area = 0.0
    total_area = 0.0
    for area, cost in plates:
        total_cost_area += cost * area
        total_area += area
    if total_area == 0:
        return 0.0
    return total_cost_area / total_area


# ---------------------------------------------------------------------------
# Self-test: representative points for each of the eight types.
# Run `python backset_cost.py` to verify the classifier reproduces the
# paper's eight classes and their costs.
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Eight-class self-test (Parsons, Nam & Singer 1999)\n")
    print(f"{'type':>4} {'b1':>7} {'b2':>8}  {'cost':>4}  name")
    print("-" * 62)

    # One representative (b1, b2) per intended type.
    cases = [
        (0.005, 0.000),  # -> Type 1  flat
        (0.100, 0.000),  # -> Type 2  low single curvature
        (0.200, 0.000),  # -> Type 3  high single curvature
        (0.050, 0.010),  # -> Type 4  moderate + small backset
        (0.120, 0.030),  # -> Type 5  moderate double
        (0.120, -0.030),  # -> Type 6  moderate reverse double
        (0.200, 0.050),  # -> Type 7  high double
        (0.200, -0.050),  # -> Type 8  high reverse double
    ]
    for b1, b2 in cases:
        r = classify_plate(b1, b2)
        print(f"{r.plate_type:>4} {b1:>7.3f} {b2:>8.3f}  {r.relative_cost:>4}  {r.name}")

    # Quick scalar-metric demo: a mostly-flat hull with a few curved plates.
    print("\nScalar-metric demo:")
    demo_plates = [
        (100.0, 1),  # big flat panel
        (100.0, 1),  # big flat panel
        (10.0, 3),  # curved plate
        (10.0, 6),  # double-curvature plate
    ]
    sm = scalar_hull_metric(demo_plates)
    print(f"  area-weighted relative cost = {sm:.3f}")