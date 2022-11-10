'''
Tests (pytest) for msh2vtu
'''
from context import msh2vtu
import os
import argparse   
import meshio
import warnings

parser = argparse.ArgumentParser()
working_dir=os.path.dirname(__file__)


def test_msh_vtu():
    '''
    Test whether msh2vtu finishes without errors 
    and generated vtu-files are readable.

    Returns
    -------
    None.

    '''
    msh_files = ["cube_tet.msh", "square_quad.msh", "square_tri.msh", "square_with_circular_hole.msh"]
    vtu_files = ["cube_tet_boundary.vtu",
                 "cube_tet_domain.vtu",
                 "cube_tet_physical_group_A.vtu",
                 "cube_tet_physical_group_B.vtu",
                 "cube_tet_physical_group_vier.vtu",
                 "cube_tet_physical_group_Wuerfel.vtu",
                 "square_quad_boundary.vtu",
                 "square_quad_domain.vtu",
                 "square_quad_physical_group_Einheitsquadrat.vtu",
                 "square_quad_physical_group_links.vtu",
                 "square_quad_physical_group_oben.vtu",
                 "square_quad_physical_group_rechts.vtu",
                 "square_quad_physical_group_unten.vtu",
                 "square_tri_boundary.vtu",
                 "square_tri_domain.vtu",
                 "square_tri_physical_group_Einheitsquadrat.vtu",
                 "square_tri_physical_group_links.vtu",
                 "square_tri_physical_group_oben.vtu",
                 "square_tri_physical_group_rechts.vtu",
                 "square_tri_physical_group_unten.vtu",
                 "square_with_circular_hole_boundary.vtu",
                 "square_with_circular_hole_domain.vtu",
                 "square_with_circular_hole_physical_group_Inner_Boundary.vtu",
                 "square_with_circular_hole_physical_group_Outer_Boundary_bottom.vtu",
                 "square_with_circular_hole_physical_group_Outer_Boundary_left.vtu",
                 "square_with_circular_hole_physical_group_Outer_Boundary_right.vtu",
                 "square_with_circular_hole_physical_group_Outer_Boundary_top.vtu",
                 "square_with_circular_hole_physical_group_SG.vtu"]
    
    for msh_file in msh_files:
        args = argparse.Namespace(filename=msh_file, output='', dim=0, delz=False, swapxy=False, rdcd=True, ogs=True, ascii=False)   
        assert msh2vtu.run(args) == 0, "msh2vtu finished with errors."
    
    for vtu_file in vtu_files:
        ErrorCode = 0
        try:
            mesh = meshio.read(vtu_file)
        except:
            ErrorCode = 1
            #print("ERROR: " + vtu_file)
        assert ErrorCode == 0, "Generated vtu-files are erroneous."
    