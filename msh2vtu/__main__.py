import argparse
from msh2vtu import run, msh2vtu_version, tested_gmsh_version, tested_meshio_version, first_meshio_version_without_meshtools

if __name__ == "__main__":
    ''' command line use '''
    
    # parsing command line arguments
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
        description = "Prepares a Gmsh-mesh for use in OGS by extracting domain-, boundary- and physical group-submeshes, and saves them in vtu-format. Note that all mesh entities should belong to physical groups.",
        epilog = "Tested with Meshio "
        + tested_meshio_version
        + " and Gmsh "
        + tested_gmsh_version
        + ". Check for changes between versions, if there are errors.",
    )
    parser.add_argument("filename", help="Gmsh mesh file (*.msh) as input data")
    parser.add_argument(
        "-g",
        "--ogs",
        action = "store_true",
        help = 'rename "gmsh:physical" to "MaterialIDs" for domains and change type of corresponding cell data to INT32',
    )
    parser.add_argument(
        "-r",
        "--rdcd",
        action = "store_true",
        help = "renumber domain cell data, physical IDs (cell data) of domains get numbered beginning with zero",
    )
    parser.add_argument(
        "-a",
        "--ascii",
        action = "store_true",	
        help = "save output files (*.vtu) in ascii format",
    )
    parser.add_argument(
        "-d",
        "--dim",
        type = int,
        default = 0,
        help = "spatial dimension (1, 2 or 3), trying automatic detection, if not given",
    )
    parser.add_argument(
        "-o",
        "--output",
        default = "",
        help = "basename of output files; if not given, then it defaults to basename of inputfile",
    )
    parser.add_argument(
        "-z",
        "--delz",
        action = "store_true",
        help = "deleting z-coordinate, for 2D-meshes with z=0, note that vtu-format requires 3D points",
    )
    parser.add_argument(
        "-s",
        "--swapxy",
        action = "store_true",
        help = "swap x and y coordinate",
    )
    parser.add_argument('-v', '--version', action='version', version='msh2vtu {} (Dominik Kern)'.format(msh2vtu_version)) 
        
    args = parser.parse_args()   
    
    run(args)
