#!/usr/bin/env python3
"""
generate_meshes.py — Generate benchmark meshes using Gmsh
Generates SENT and TPB meshes at multiple resolutions.

Requirements: pip install gmsh
Paper Appendix D mesh specifications:
  - SENT: L×H = 1.0×0.5 m, h_bulk=0.010 m, notch at left edge, a0=0.2 m
  - TPB:  L×H = 4.0×1.0 m, h_bulk=0.020 m, notch at bottom center, a0=0.3 m
"""

import sys
import os
try:
    import gmsh
except ImportError:
    print("ERROR: gmsh not installed. Run: pip install gmsh")
    sys.exit(1)


def generate_sent(h_bulk=0.010, output_path="meshes/sent_h010.msh"):
    """Single Edge Notched Tension specimen."""
    gmsh.initialize()
    gmsh.model.add("SENT")
    gmsh.option.setNumber("General.Terminal", 0)  # suppress output

    L = 1.0   # width [m]
    H = 0.5   # height [m]
    a0 = 0.2  # notch length [m]

    # Domain boundary
    p1 = gmsh.model.geo.addPoint(0,   0,  0, h_bulk)
    p2 = gmsh.model.geo.addPoint(L,   0,  0, h_bulk)
    p3 = gmsh.model.geo.addPoint(L,   H,  0, h_bulk)
    p4 = gmsh.model.geo.addPoint(0,   H,  0, h_bulk)

    # Notch tip
    pn = gmsh.model.geo.addPoint(a0,  H/2, 0, h_bulk*0.3)  # refined at tip

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, pn)   # top-left to notch tip
    l5 = gmsh.model.geo.addLine(pn, p1)   # notch tip to bottom-left

    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.synchronize()

    # Mesh size at notch tip (h_tip = 0.3*h_bulk per paper)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", [pn])

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", h_bulk*0.3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", h_bulk)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 3*h_bulk)

    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.Algorithm", 5)  # Delaunay

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(1)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    gmsh.write(output_path)
    n = len(gmsh.model.mesh.getNodes()[0])
    e = len([t for t in gmsh.model.mesh.getElements(2)[1][0]])
    gmsh.finalize()
    print(f"  SENT (h={h_bulk}): {n} nodes, {e} elements → {output_path}")
    return output_path


def generate_tpb(h_bulk=0.020, output_path="meshes/tpb_h020.msh"):
    """Three-Point Bending specimen."""
    gmsh.initialize()
    gmsh.model.add("TPB")
    gmsh.option.setNumber("General.Terminal", 0)

    L = 4.0   # length [m]
    H = 1.0   # height [m]
    a0 = 0.3  # notch depth [m]

    # Domain
    p1 = gmsh.model.geo.addPoint(0,   0,   0, h_bulk)
    p2 = gmsh.model.geo.addPoint(L,   0,   0, h_bulk)
    p3 = gmsh.model.geo.addPoint(L,   H,   0, h_bulk)
    p4 = gmsh.model.geo.addPoint(0,   H,   0, h_bulk)

    # Notch at bottom center
    xmid = L/2
    pn = gmsh.model.geo.addPoint(xmid, a0, 0, h_bulk*0.3)  # notch tip
    pb1 = gmsh.model.geo.addPoint(xmid-0.001, 0, 0, h_bulk*0.3)
    pb2 = gmsh.model.geo.addPoint(xmid+0.001, 0, 0, h_bulk*0.3)

    l1  = gmsh.model.geo.addLine(p1, pb1)
    ln1 = gmsh.model.geo.addLine(pb1, pn)
    ln2 = gmsh.model.geo.addLine(pn, pb2)
    l2  = gmsh.model.geo.addLine(pb2, p2)
    l3  = gmsh.model.geo.addLine(p2, p3)
    l4  = gmsh.model.geo.addLine(p3, p4)
    l5  = gmsh.model.geo.addLine(p4, p1)

    loop = gmsh.model.geo.addCurveLoop([l1, ln1, ln2, l2, l3, l4, l5])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.synchronize()

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", [pn])

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", h_bulk*0.3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", h_bulk)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 3*h_bulk)

    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.Algorithm", 5)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(1)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    gmsh.write(output_path)
    gmsh.finalize()
    print(f"  TPB (h={h_bulk}): → {output_path}")
    return output_path


def main():
    print("Generating benchmark meshes...")
    os.makedirs("meshes", exist_ok=True)

    # SENT meshes at 4 refinement levels (for convergence study)
    for h in [0.040, 0.020, 0.010, 0.005]:
        tag = str(h).replace(".", "")
        generate_sent(h_bulk=h, output_path=f"meshes/sent_h{tag}.msh")

    # TPB meshes
    for h in [0.040, 0.020, 0.010]:
        tag = str(h).replace(".", "")
        generate_tpb(h_bulk=h, output_path=f"meshes/tpb_h{tag}.msh")

    print("Done. Meshes written to meshes/")


if __name__ == "__main__":
    main()
