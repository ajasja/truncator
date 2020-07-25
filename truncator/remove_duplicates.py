from .pymol_utils import *
import os

def remove_duplicates_rmsd(pdb_list, rmsd_cutoff=0.5, extra_list=None, verbose=False):
    """Remove duplicates from a list of pdbs based on the RMSD. Eveythin below the rmsd_cutoff is considert a duplicate. 
    The first file name is kept. Returns a list of unique structures. Takes an extra list (of indices etc, than gets trited the same way)"""
    m=0
    
    pdb_list = list(pdb_list)
    if extra_list:
        assert len(extra_list)==len(pdb_list), 'Lengths of pdb_list and extra list must match'
        extra_list=list(extra_list)
    if len(pdb_list)<2:
        if extra_list:
            return pdb_list, extra_list
        return pdb_list 

    while m<len(pdb_list):
        n=m+1
        while n<len(pdb_list):
            if verbose:
                print(m,n)
            assert os.path.exists(pdb_list[m]), f'{pdb_list[m]} does not exist' 
            assert os.path.exists(pdb_list[n]), f'{pdb_list[n]} does not exist'
            rmsd = get_rmsd(pdb_list[m], pdb_list[n], align=False, cleanup=True)
            if verbose:
               print(f"{pdb_list[m]}\n{pdb_list[n]}\nRMDS:{rmsd}") 
            if rmsd<rmsd_cutoff:
                if verbose:
                    print(f'Deleting index {n}/{len(pdb_list)}, {pdb_list[n]}') 
                del pdb_list[n]
                if extra_list:
                    del extra_list[n]
                
                n=n-1 #Adjust index, since all are moved up. #IMPORTAN BUG FIX, without this it's skipping some structures

            n=n+1
        m=m+1
    
    if extra_list:
        return pdb_list, extra_list
    return pdb_list

