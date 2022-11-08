'''
Running msh2vtu on standard meshes (msh-files) and then checking if generated
vtu-files can be read.
'''
from context import msh2vtu
import os
import argparse   
import meshio

parser = argparse.ArgumentParser()
working_dir=os.path.dirname(__file__)


def test_msh():
    msh_files = [x for x in os.listdir(working_dir) if x.endswith(".msh")]
    for msh_file in msh_files:
        args = argparse.Namespace(filename=msh_file, output='', dim=0, delz=False, swapxy=False, rdcd=True, ogs=True, ascii=False)   
        assert msh2vtu.run(args) == 0
    

def test_vtu():
    vtu_files = [x for x in os.listdir(working_dir) if x.endswith(".vtu")]
    for vtu_file in vtu_files:
        ErrorCode = 0
        try:
            mesh=meshio.read(vtu_file)
        except:
            ErrorCode = 1
            #print("ERROR: " + vtu_file)
        assert ErrorCode == 0