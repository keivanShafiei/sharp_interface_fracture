#!/usr/bin/env python3
"""
Generate TPB (Three-Point Bending) mesh using Gmsh
Compatible with sharp_interface_fracture
"""

import gmsh
import sys

def create_tpb_mesh(output_file="tpb_rect_coarse2.msh", mesh_size=0.05):
    """
    Create a TPB geometry:
    - Rectangle: L=4.0, H=1.0 (typical TPB aspect ratio)
    - Notch at bottom: depth=0.5H
    - Support points at bottom corners
    - Load point at top center
    """
    
    gmsh.initialize()
    gmsh.model.add("TPB")
    
    # Geometry parameters
    L = 4.0      # Length
    H = 1.0      # Height
    notch_depth = 0.5 * H  # Notch depth (50% of height)
    
    lc = mesh_size  # Characteristic length (controls mesh size)
    
    # Points for outer rectangle
    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc)
    p2 = gmsh.model.geo.addPoint(L, 0.0, 0.0, lc)
    p3 = gmsh.model.geo.addPoint(L, H, 0.0, lc)
    p4 = gmsh.model.geo.addPoint(0.0, H, 0.0, lc)
    
    # Notch points (at bottom center)
    notch_x = L / 2.0
    p5 = gmsh.model.geo.addPoint(notch_x - 0.05, H - notch_depth, 0.0, lc * 0.5)  # Notch tip
    p6 = gmsh.model.geo.addPoint(notch_x + 0.05, H - notch_depth, 0.0, lc * 0.5)  # Notch tip
    p7 = gmsh.model.geo.addPoint(notch_x - 0.1, H, 0.0, lc)  # Notch left
    p8 = gmsh.model.geo.addPoint(notch_x + 0.1, H, 0.0, lc)  # Notch right
    
    # Lines for outer boundary
    l1 = gmsh.model.geo.addLine(p1, p2)  # Bottom
    l2 = gmsh.model.geo.addLine(p2, p3)  # Right
    l3_top = gmsh.model.geo.addLine(p3, p8)  # Top right (before notch)
    l3_notch_right = gmsh.model.geo.addLine(p8, p6)  # Notch right slope
    l3_notch_bottom = gmsh.model.geo.addLine(p6, p5)  # Notch bottom
    l3_notch_left = gmsh.model.geo.addLine(p5, p7)  # Notch left slope
    l3_top_left = gmsh.model.geo.addLine(p7, p4)  # Top left (after notch)
    l4 = gmsh.model.geo.addLine(p4, p1)  # Left
    
    # Create surface
    loop1 = gmsh.model.geo.addCurveLoop([
        l1, l2, l3_top, l3_notch_right, l3_notch_bottom, l3_notch_left, l3_top_left, l4
    ])
    surf = gmsh.model.geo.addPlaneSurface([loop1])
    
    # Synchronize geometry
    gmsh.model.geo.synchronize()
    
    # Add physical tags for boundary conditions
    # Support points (bottom)
    gmsh.model.addPhysicalGroup(0, [p1], 1, name="Support_Left")
    gmsh.model.addPhysicalGroup(0, [p2], 2, name="Support_Right")
    
    # Load point (top center)
    gmsh.model.addPhysicalGroup(0, [p4], 3, name="Load_Point")
    
    # Crack surfaces (notch region)
    gmsh.model.addPhysicalGroup(1, [l3_notch_left, l3_notch_bottom, l3_notch_right], 4, name="Initial_Crack")
    
    # Material region
    gmsh.model.addPhysicalGroup(2, [surf], 5, name="Material")
    
    # Generate 2D mesh
    gmsh.model.mesh.generate(2)
    
    # Optimize mesh
    gmsh.model.mesh.optimize("Netgen")
    
    # Write mesh file
    gmsh.write(output_file)
    print(f"✓ Mesh created: {output_file}")
    
    # Get mesh statistics
    gmsh.model.mesh.getStatistics()
    
    gmsh.finalize()
    return output_file

if __name__ == "__main__":
    mesh_size = 0.05 if len(sys.argv) < 2 else float(sys.argv[1])
    output = "tpb_rect_coarse2.msh" if len(sys.argv) < 3 else sys.argv[2]
    
    print(f"Creating TPB mesh with element size ≈ {mesh_size}")
    create_tpb_mesh(output, mesh_size)
    print(f"✓ Done! Mesh saved to: {output}")
