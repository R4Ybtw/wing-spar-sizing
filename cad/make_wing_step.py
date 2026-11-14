"""Build the half-wing (for Flow Simulation) and the spar tube (for SolidWorks
Simulation) from data/aircraft.json, both split into spanwise strips.

    pip install cadquery
    python make_wing_step.py                  default tube 25 x 0.7 mm
    python make_wing_step.py 25 1.0           any tube: outer diameter, wall (mm)

Coordinates (mm): span along +Z from the root (z = 0, the symmetry plane) to
the tip; chord along +X with the QUARTER CHORD on the Z axis; lift along +Y.
Each strip is a separate piece fused WITHOUT merging faces, so every strip
keeps its own faces: one Flow Simulation goal per strip, one Simulation load
per strip. If SolidWorks merges them on import, re-split with Split Line
(planes at the z values the script prints).
"""
import json
import sys
import cadquery as cq

with open("../data/aircraft.json") as f:
    ac = json.load(f)

c = ac["chord"] * 1000
s = ac["span"] / 2 * 1000
n = ac["n_segments"]
code = ac["airfoil"]
m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
D = float(sys.argv[1]) if len(sys.argv) > 2 else 25.0
tw = float(sys.argv[2]) if len(sys.argv) > 2 else 0.7


def naca4(npts=81):
    """Upper and lower surface points, LE -> TE, closed sharp TE (0.1036 term)."""
    import math
    up, lo = [], []
    for i in range(npts):
        x = 0.5 * (1 - math.cos(math.pi * i / (npts - 1)))
        yt = 5 * t * (0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
        if m > 0 and x < p:
            yc, dy = m / p**2 * (2 * p * x - x**2), 2 * m / p**2 * (p - x)
        elif m > 0:
            yc, dy = m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x**2), 2 * m / (1 - p) ** 2 * (p - x)
        else:
            yc, dy = 0.0, 0.0
        th = math.atan(dy)
        up.append(((x - yt * math.sin(th) - 0.25) * c, (yc + yt * math.cos(th)) * c))
        lo.append(((x + yt * math.sin(th) - 0.25) * c, (yc - yt * math.cos(th)) * c))
    lo[-1] = up[-1]          # exactly closed trailing edge
    return up, lo


up, lo = naca4()
edges = [s * k / n for k in range(n + 1)]


def section_face(z):
    upper = cq.Edge.makeSpline([cq.Vector(x, y, z) for x, y in reversed(up)])
    lower = cq.Edge.makeSpline([cq.Vector(x, y, z) for x, y in lo])
    return cq.Face.makeFromWires(cq.Wire.assembleEdges([upper, lower]))


def strips(make_piece):
    pieces = [make_piece(edges[k], edges[k + 1]) for k in range(n)]
    body = pieces[0]
    for pc in pieces[1:]:
        body = body.fuse(pc, glue=True)      # glue: keep the faces at each strip boundary
    return body


wing = strips(lambda z0, z1: cq.Solid.extrudeLinear(section_face(z0), cq.Vector(0, 0, z1 - z0)))
tube = strips(lambda z0, z1: cq.Solid.makeCylinder(D / 2, z1 - z0, cq.Vector(0, 0, z0))
              .cut(cq.Solid.makeCylinder(D / 2 - tw, z1 - z0, cq.Vector(0, 0, z0))))

cq.exporters.export(cq.Workplane().add(wing), "half_wing_strips.step")
cq.exporters.export(cq.Workplane().add(tube), f"spar_tube_{D:g}x{tw:g}_strips.step")
for name, shp in (("half wing", wing), ("spar tube", tube)):
    bb = shp.BoundingBox()
    print(f"{name}: {bb.xlen:.1f} x {bb.ylen:.1f} x {bb.zlen:.1f} mm, faces {len(shp.Faces())}, "
          f"volume {shp.Volume()/1e3:.1f} cm^3, valid={shp.isValid()}")
print("strip boundaries z [mm]:", ", ".join(f"{e:.0f}" for e in edges))
