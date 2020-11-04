#!/usr/bin/env python
# TODO extension to 3D, also option to save 2D as [x,y] without z=0

import meshio
import os
import numpy
import argparse
import warnings


def line_mesh_prune(points, input_cells):  # remove orphaned points
    # "old" means from the input mesh and "new" the mesh of connected points only
    original_shape = input_cells.shape  # for reconstrution after flatten
    old_points = input_cells.flatten()  # 1d-array needed
    unique_points, unique_inverse = numpy.unique(old_points, return_inverse=True)
    new_points = points[unique_points]  # extract only used nodes
    new_cells = unique_inverse.reshape(original_shape)  # update cell connectivity
    return new_points, new_cells


# parsing command line arguments
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(
    description="Prepare Gmsh-mesh for use in OGS by extracting domain-, boundary- and physical group-meshes and save them in vtu-format. Note that all mesh entities must belong to a physical group!",
    epilog="No cell data are written for boundaries (lines).",
)
parser.add_argument("filename", help="Gmsh mesh file (*.msh) as input data.")
parser.add_argument(
    "--renumber",
    action="store_true",
    help="Renumber physical IDs of domains starting by zero (boundary IDs are ignored).",
)
parser.add_argument(
    "--ogs",
    action="store_true",
    help='Rename "gmsh:physical" to "MaterialIDs" for domains and change type of corresponding cell data to INT32.',
)
parser.add_argument(
    "-a",
    "--ascii",
    action="store_true",
    help="Save output files (*.vtu) in ascii format.",
)
parser.add_argument(
    "-o",
    "--output",
    default="",
    help="Base name of output files; if not given, then it defaults to basename of inputfile.",
)

args = parser.parse_args()


# check if input file exists and is in gmsh-format
if os.path.isfile(args.filename):
    filename_without_extension = os.path.splitext(args.filename)[0]
    file_extension = os.path.splitext(args.filename)[1]
    if file_extension != ".msh":
        warnings.warn("Warning, input file seems not to be in gmsh-format (*.msh)")
else:
    raise FileNotFoundError


# derive output filenames
if args.output == "":  # no parameter given, use same basename as input file
    output_basename = filename_without_extension
else:
    output_basename = args.output


# read in mesh
mesh = meshio.read(args.filename)
# points: 2Darray( point, xyz)
# cells: list( cellblocks( type='line'/'triangle', data=2Darray(element, points) ))
# cell_data: dict('gmsh:physical', 2Darray(geo_tag, ph_tag))
# field_data: dict('top'/'bottom'/PH_NAME, 1Darray (ph_tag, geo_type)
points, cells_dict = mesh.points, mesh.cells_dict
cell_data_dict, field_data = mesh.cell_data_dict, mesh.field_data


# some variable declarations
ph_index = 0  # to access physical id in field data
geo_index = 1  # to access geometrical id in field data
line_id = 1  # geometry type
triangle_id = 2  # geometry type
gmshdict = {line_id: "line", triangle_id: "triangle"}  # gmsh convention
gmsh_string = "gmsh:physical"
ogs_string = "MaterialIDs"

# if user wants to change 'gmsh:physical" to "MaterialIDs" or not
if args.ogs:
    # boundary_data_string=ogs_string # do not write MaterialID for boundaries
    domain_data_string = ogs_string
    selection_data_string = ogs_string  # only used for domains
else:
    # boundary_data_string=gmsh_string
    domain_data_string = gmsh_string
    selection_data_string = gmsh_string  # only used for domains

# if user wants physical group numbering of domains beginning with zero
id_offset = 0  # initial value, zero will not change anything
if args.renumber:
    # find minimum physical_id of domains (triangles)
    id_list_domains = []
    for dataset in field_data.values():  # go through all physical groups
        if (
            dataset[geo_index] == triangle_id
        ):  # only for domains (triangles), ignore boundaries (lines)
            id_list_domains.append(dataset[0])  # append physical id
    if len(id_list_domains):  # if there are some domains..
        id_offset = min(id_list_domains)  # ..then find minimal physical id


# extract domain and boundary mesh

domain_cells = mesh.cells_dict[gmshdict[triangle_id]]
domain_cell_data = mesh.cell_data_dict[gmsh_string][gmshdict[triangle_id]]
if args.ogs:
    domain_cell_data = numpy.int32(domain_cell_data)  # ogs needs MaterialIDs as int32

# write results to file
if len(domain_cells) and len(domain_cells) == len(
    domain_cell_data
):  # only if there are (consistent) data to write
    domain_mesh = meshio.Mesh(
        points=points,
        cells=[(gmshdict[triangle_id], domain_cells)],
        cell_data={domain_data_string: [domain_cell_data - id_offset]},
    )
    domain_mesh.prune()
    meshio.write(output_basename + "_domain.vtu", domain_mesh, binary=not args.ascii)
else:
    print("No or inconsistent domain-cells found, no domain-mesh written.")

boundary_cells = mesh.cells_dict[gmshdict[line_id]]
# write results to file
if len(boundary_cells):  # only if there are data to write
    boundary_mesh = meshio.Mesh(
        points=points, cells=[(gmshdict[line_id], boundary_cells)]
    )
    # boundary_mesh.prune()   # strange, leads to error
    meshio.write(
        output_basename + "_boundary.vtu", boundary_mesh, binary=not args.ascii
    )
else:
    print("No boundary-cells found, no boundary-mesh written.")

# Now we want to extract subdomains given by physical groups in gmsh,
# so we need an additional loop.

# name=user-defined name of physical group, data=[physical_id, geometry_id]
for name, data in field_data.items():

    ph_id = data[ph_index]  # selection by physical id (user defined)
    geo_id = data[geo_index]  # 1 or 2
    cell_type = gmshdict[geo_id]  # 'line' or 'triangle'
    selection_index = cell_data_dict[gmsh_string][cell_type] == ph_id
    selected_cells = cells_dict[cell_type][selection_index]
    selected_cell_data = cell_data_dict[gmsh_string][cell_type][selection_index]
    if args.ogs:
        selected_cell_data = numpy.int32(selected_cell_data)

    if len(selected_cells):  # if there are some data
        if data[geo_index] == line_id:  # prune line mesh and create submesh
            new_points, new_cells = line_mesh_prune(points, selected_cells)
            submesh = meshio.Mesh(points=new_points, cells=[(cell_type, new_cells)])
        elif (
            data[geo_index] == triangle_id
        ):  # create submesh, including cell_data, and prune it
            submesh = meshio.Mesh(
                points=points,
                cells=[(cell_type, selected_cells)],
                cell_data={selection_data_string: [selected_cell_data - id_offset]},
            )
            submesh.prune()
        else:
            print("Unknown geometry id encountered, empty submesh.")
            submesh = meshio.Mesh(points=[], cells=[])
        outputfilename = output_basename + "_physical_group_" + name + ".vtu"
        meshio.write(outputfilename, submesh, binary=not args.ascii)
    else:
        print("No cells found for physical group " + name + ", no submesh written.")
