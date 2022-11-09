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


def test_msh():
    '''
    Test whether msh2vtu finishes without errors.

    Returns
    -------
    None.

    '''
    vtu_files = [x for x in os.listdir(working_dir) if x.endswith(".vtu")]
    if len(vtu_files) > 0:
        warnings.warn("Warning, there are some vtu-files in this folder before msh2vtu was run.", stacklevel=2)
    
    msh_files = [x for x in os.listdir(working_dir) if x.endswith(".msh")]
    assert len(msh_files) > 0, "There are no msh-files to test."
    
    for msh_file in msh_files:
        args = argparse.Namespace(filename=msh_file, output='', dim=0, delz=False, swapxy=False, rdcd=True, ogs=True, ascii=False)   
        assert msh2vtu.run(args) == 0, "msh2vtu finished with errors."
    

def test_vtu():
    '''
    Test whether generated files can be read.

    Returns
    -------
    None.

    '''
    vtu_files = [x for x in os.listdir(working_dir) if x.endswith(".vtu")]
    assert len(vtu_files) > 0, "There are no results to test."
    
    for vtu_file in vtu_files:
        ErrorCode = 0
        try:
            mesh = meshio.read(vtu_file)
        except:
            ErrorCode = 1
            #print("ERROR: " + vtu_file)
        assert ErrorCode == 0, "Generated vtu-files are erroneous."