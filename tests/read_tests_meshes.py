import os
import sys
import meshio

path="/home/dominik/git_repos/msh2vtu/tests/"
vtu_files = [x for x in os.listdir(path) if x.endswith(".vtu")]

error_code=0
for vtu_file in vtu_files:
    try:
        mesh=meshio.read(path+vtu_file)
    except:
        error_code=1
        print("ERROR: " + vtu_file)

if error_code==1:
    sys.exit(1)
else:
    sys.exit(0)

