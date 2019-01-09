import truncator
from .utils import find_helices, python_to_pymol_range, list_to_pymol_sel_str, flatten_list
import itertools

import numpy
import numpy as np

from glob import glob
import os

import pymol
from pymol import cmd


def cut_bundles(struct_name, out_dir, tol_A=0.5, num_heptads=3, step_heptad_fraction=1/4*0.9, cmd=None):
    """Cuts bundles to a number of heptades in steps of heptad fractions"""
    md = {}
    md['cutting.tol_A']=tol_A
    md['cutting.num_heptads']=num_heptads
    md['cutting.step_heptad_fraction']=step_heptad_fraction
    md['cutting.base_struct_name']=struct_name
    md['cutting.base_struct_name_full']=os.path.abspath(struct_name)

    if cmd is None:
        import pymol
        cmd =  pymol.cmd
        cmd.do("delete all")
    base_name = truncator.basename_noext(struct_name)
    cmd.load(struct_name, object=base_name)
    model = cmd.get_model("name ca")
    sec_struct = [at.ss for at in model.atom]
    
    #find helix ends
    helices_pos = find_helices(sec_struct)
    helices_pos_pymol = [python_to_pymol_range(*pos) for pos in helices_pos]
    helices_ends_pymol = list_to_pymol_sel_str(flatten_list(helices_pos_pymol))
    
    #rechain helices
    chainIDs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    for n, (_from, _to) in enumerate(helices_pos_pymol):
        #print(n,_from, _to )
        cmd.alter(f"resi {_from}-{_to}", f"chain='{chainIDs[n]}'")
        
    #get helix orientations as a dictionary 
    helix_orientations = {}
    for n, (_from, _to) in enumerate(helices_pos_pymol):
        sel_str =  f"resi {_from}+{_to}"
        #print(n,_from, _to, sel_str)
        
        helix_coords = cmd.get_coords(sel_str)

        z_diff = np.diff(np.array(helix_coords)[:,2])[0]
        if z_diff > 0:
            helix_orientations[chainIDs[n]] = "U"
        else:
            helix_orientations[chainIDs[n]] = "D"
            
    #find bottom parts
    lower_coords = cmd.get_coords(f"name ca and z<0 and resi {helices_ends_pymol}")
    upper_coords = cmd.get_coords(f"name ca and z>0 and resi {helices_ends_pymol}")
    #get maximum (highest point) on lower end
    highest_bottom = np.max(lower_coords[:,2])
    #and lowest point on top
    lowest_top = np.min(upper_coords[:,2])
    
    #measure the 7 residue rise 
    _start=helices_pos_pymol[0][0]
    _end = _start+7-1
    heptad_delta_z = np.diff(np.array(cmd.get_extent(f"name ca and resi {_start}-{_end}"))[:,2])[0]
    
    step_delta_z=heptad_delta_z*step_heptad_fraction
    height = heptad_delta_z * num_heptads
    
    md['cutting.heptad_delta_z'] = heptad_delta_z
    md['cutting.helix_orientations'] = helix_orientations
    
    bottoms = np.arange(highest_bottom, lowest_top, step_delta_z)
    #filter the bottoms that would not result in a full length protein
    bottoms = [bottom for bottom in bottoms if (bottom+height)< (lowest_top+tol_A)]
    
    #output the files
    truncator.make_dirs(out_dir)
    files = []
    for bottom in bottoms:
        _from =  '%05.2f' % (bottom - tol_A)
        _to   =  '%05.2f' % (bottom + height + tol_A)
        out_name = f"{out_dir}/{base_name}__numH{num_heptads}__from{_from}__to{_to}.pdb"
        md['cutting.from']=_from
        md['cutting.to']=_to
        
        sel_str = f"byres (name ca and ss H and z>{_from} and z<{_to})"
        print(out_name, ": ", sel_str)
        cmd.save(out_name, sel_str)
        truncator.write_json(truncator.replace_extension(out_name,'.info') , md)
        files.append(out_name)

    return files


def regroup_chains(struct_name, out_dir, new_chain_A, new_chain_B, cmd=None):
    """Given a PDB with multiple chains, it changes the chains to A and B for interface evaluation. 
    Takes a PDB and outputs a PDB.
    """
    if cmd is None:
        import pymol
        cmd =  pymol.cmd
        cmd.do("delete all")

    base_name = truncator.basename_noext(struct_name)
    md = truncator.read_info_file(struct_name,".info")
    md['regroup.new_chain_A'] = new_chain_A
    md['regroup.new_chain_B'] = new_chain_B
    cmd.do("delete all")
    cmd.load(struct_name, object=base_name)

    #A
    for old_chain in new_chain_A:
        cmd.alter(f"chain {old_chain}", f"chain='A'")

    #B
    for old_chain in new_chain_B:
        cmd.alter(f"chain {old_chain}", f"chain='B'")

    truncator.make_dirs(out_dir)
    out_name = f"{out_dir}/{base_name}__gr{new_chain_A}-{new_chain_B}.pdb"

    print(out_name)
    cmd.save(out_name)
    truncator.write_json(out_name.replace(".pdb",".info"), md)
    return out_name
