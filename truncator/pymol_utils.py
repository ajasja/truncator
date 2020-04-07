try:
    import pymol
    from pymol import cmd
    from pymol import stored
    from pymol import CmdException
except:
    print("Pymol is most likely not installed")


def get_ss(sel_str):
    """returns the secondary structure of selection as a string"""
    stored.ss_array = []
    cmd.iterate(f"{sel_str} and (name CA)", "stored.ss_array.append(ss)")
    return "".join(stored.ss_array)


def get_selection_property(sel_str, property="resi"):
    """returns the secondary structure of selection as a string"""
    stored.seli_array = []
    cmd.iterate(f"{sel_str} and (name CA)", f"stored.seli_array.append({property})")
    return stored.seli_array


get_sp = get_selection_property


import itertools



def get_helices(sele="all"):
    helices = list()
    prevss = "nope"
    for atom in cmd.get_model(sele + " and name CA").atom:
        # print(atom.ss, atom.resi, atom.resi_number, atom.name)
        if atom.ss == "H":
            if atom.ss != prevss:
                helices.append(list())
            helices[-1].append(atom.resi)
        prevss = atom.ss
    return helices


def sel_first_helix(n, sele="all", onlyh=True, sele_name="sele"):
    helices = get_helices(sele)
    if int(onlyh):
        hresi = itertools.chain(*helices[: int(n)])
        sele_new = "resi " + "+".join(r for r in hresi)
    else:
        sele_new = "resi 0-%s" % helices[int(n) - 1][-1]

    sele_new = sele + " AND (" + sele_new + ")"
    cmd.select(sele_name, selection=sele_new)
    return sele_new


def sel_last_helix(n, sele="all", onlyh=True, sele_name="sele"):
    helices = get_helices(sele)
    if int(onlyh):
        hresi = itertools.chain(*helices[-int(n) :])
        sele_new = "resi " + "+".join(r for r in hresi)
    else:
        sele_new = "resi %s-99999" % helices[-int(n)][0]

    sele_new = sele + " AND (" + sele_new + ")"
    cmd.select(sele_name, selection=sele_new)
    return sele_new


def sel_terminal_helix(n, sele="all", onlyh=True, sele_name="sele"):
    if int(n) > 0:
        return sel_first_helix(int(n), sele=sele, onlyh=onlyh, sele_name=sele_name)
    else:
        return sel_last_helix(-int(n), sele=sele, onlyh=onlyh, sele_name=sele_name)


cmd.extend("sel_first_helix", sel_first_helix)
cmd.extend("sel_last_helix", sel_last_helix)
cmd.extend("sel_terminal_helix", sel_terminal_helix)


def clash_check_CA(selA, selB, distance=4.0, sele_name=None):
    """Returns the number of CA atoms in selA that lie closer than radii to any CA atom in selB"""
    ##
    selA = f"({selA}) and name CA"
    selB = f"({selB}) and name CA"

    sel_str = f"(({selA}) within {distance} of ({selB}))"
    model = cmd.get_model(sel_str)
    if not sele_name is None:
        cmd.select(sele_name, selection=sel_str)
    return model.nAtom


def get_distance(c1, c2):
    return ((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2) ** 0.5


def get_alignment_map(from_sel, to_sel, max_distance=1.8):
    """Finds atoms in to_sel that are close to from_sel"""
    from_sel = f"({from_sel}) and name CA"
    to_sel = f"({to_sel}) and name CA"
    mapping = {}
    distances = {}
    from_model = cmd.get_model(from_sel)

    for at in from_model.atom:
        to_model = cmd.get_model(
            f"({to_sel}) within {max_distance} of ({from_sel} and chain {at.chain} and resi {at.resi})"
        )
        # print(f"({to_sel}) within {max_distance} of ({from_sel} and chain {at.chain} and resi {at.resi})")
        # print(to_model.nAtom)
        if to_model.nAtom > 1:
            print(
                f"WARNING: more than one atom ({to_model.nAtom}) within {from_sel} and chain {at.chain} and resi {at.resi}"
            )
        ch_res_id = f"{at.chain}_{at.resi}"
        if to_model.nAtom > 0:
            to_model_at = to_model.atom[0]
            mapping[ch_res_id] = f"{to_model_at.chain}_{to_model_at.resi}"
            distances[ch_res_id] = get_distance(at.coord, to_model_at.coord)
        else:
            mapping[ch_res_id] = None
            distances[ch_res_id] = None
    return mapping, distances


def color_by_selector_array(name, selector, color="red"):
    if isinstance(selector, str):
        selector = selector.split(",")
    for n, i in enumerate(selector, start=1):
        if int(i) > 0:
            cmd.color(color, f"{name} and resi {n}")


def get_pymol_selector_from_rosetta_selector_str(sel_str):
    """Given a list of pdb_ids 10_A, or 10A print a pymol selector"""
    ids = sel_str.split(",")
    result = []
    # TODO group by chains?
    for id_ in ids:
        chain = id_[-1:]
        resi = id_[:-1]
        result.append(f"(resi {resi} and chain {chain})")
    result = "(" + " or ".join(result) + ")"
    return result
    
# get_pymol_selector_from_rosetta_selector_str('158B,160B,161B,162B,164B,165B,167B,168B,169B,171B,172B,194B,195B,196B,197B,198B,199B,200B,201B,202B,203B,204B,205B,206B,208B,209B')


def color_by_pdb_id_list(name, pdb_ids, color="red"):
    """Given a list of pdb_ids 10_A, or 10A print a pymol selector"""
    pass


cmd.extend("color_by_selar", color_by_selector_array)


import re


def grep_file(file_name, expression):
    """Searches file and returns lines matching regular expression"""
    result = []
    
    with open(file_name, 'r') as file:
        for line in file:
            #print(line)
            if re.search(expression, line): 
                result.append(line.strip())
    return result

#grep_file('out/repacked/ALAF05_8x/ALAF05_8x.pdb', r"^.*_pymol_selection")

def apply_read_selections(sel_str, print_cmd=False):
    """Applies read selection to pymol"""
    name, cmd_str = sel_str.split(' ', maxsplit=1)
    name = name.replace('_pymol_selection', "")
    cmd_str = cmd_str.replace('rosetta_sele', name)
    if print_cmd:
        print(cmd_str)
    else:
        cmd.do(cmd_str)
        
def import_selections_from_pdb(pdb_path, print_cmd=False):
    """Loads selection strings from the PDB file"""
    seles = grep_file(pdb_path, r"^.*_pymol_selection")
    for sel in seles:
        apply_read_selections(sel, print_cmd)


def load_with_selections(pdb_path, print_cmd=False):
    """Loads a PDB and it's selections"""
    cmd_str=f"load {pdb_path}"
    if print_cmd:
        print(cmd_str)
    else:
        cmd.do(cmd_str)
    import_selections_from_pdb(pdb_path, print_cmd)

load_ws = load_with_selections
cmd.extend("load_ws", load_with_selections)

import gzip
def read_file(filename):
    """Reads a file. Supports gzip files if ending is .gz"""
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as myfile:
            data = myfile.read()
    else:
        with open(filename, 'rt') as myfile:
            data = myfile.read()
    return data  

def apply_unsat_group(unsat_groups, print_cmd=False):
    """Parses an unsat group from BuriedUnsatHbonds (name: \n Unsatisfied HEAVY polar atom at residue 51: HIS  ND1 \n ... )"""
    lines = unsat_groups.split('\n')
    #name is the first line, skip the last colon
    name = lines.pop(0)[:-2]
    #print(name)
    selections = []
    for line in lines:
        line = line.strip()
        if line == '': continue #skip empty lines
        #split by spaces and take the last three fields
        spl = line.split()
        resnum = spl[-3][:-1] #skip the colum
        resname = spl[-2]
        atom   = spl[-1]
        
        sele = f"{resname}`{resnum}/{atom}"
        #print(sele)
        selections.append(sele)
    selections_str = " or ".join(selections)
    cmd_str = f"select {name}_unsats, {selections_str}"
    
    if print_cmd:
        print(cmd_str)
    else:
        cmd.do(cmd_str)

import re
find_unsat_sections = re.compile(r"^BuriedUnsatHbonds\s(.*?)\s\sall_heavy_atom_unsats", re.MULTILINE | re.DOTALL | re.S)

def import_unsats_from_pdb(pdb_path, print_cmd=False):
    """Greps the PDB file for unsat blocks and turns them into pymol selections"""
    file_str = read_file(pdb_path)
    unsat_groups = re.findall(find_unsat_sections, file_str)
    
    for unsat_group in unsat_groups:
        apply_unsat_group(unsat_group, print_cmd)


load_unsats = import_unsats_from_pdb
cmd.extend("load_unsats", load_unsats)

def apply_metric_from_npz(npz_file, metric_type='lddt', sel_str='all', print_cmd=False):
    '''Imports an array of values from the npz file
    Metric types are:
    lddt, pe40, pe20, pe10, pe05
    '''
    import numpy as np
    dat = np.load(npz_file)
    
    lddt = dat["lddt"] 
    esto = dat["estogram"]

    if  metric_type=='lddt':
        res = lddt
    elif metric_type=='pe40':
        res = 1-np.mean(np.sum(esto[:,:,:4], axis=-1) + np.sum(esto[:,:,11:], axis=-1), axis=-1)
    elif metric_type=='pe20':
        res = 1-np.mean(np.sum(esto[:,:,:5], axis=-1) + np.sum(esto[:,:,10:], axis=-1), axis=-1)
    elif metric_type=='pe10':
        res = 1-np.mean(np.sum(esto[:,:,:6], axis=-1) + np.sum(esto[:,:,9:], axis=-1), axis=-1)
    elif metric_type=='pe05':
        res = 1-np.mean(np.sum(esto[:,:,:7], axis=-1) + np.sum(esto[:,:,8:], axis=-1), axis=-1)
    else:
        raise f"Metric should be one of lddt, pe40, pe20, pe10, pe05, but is {metric_type}"
    
    model = cmd.get_model(f"{sel_str} and name CA")
    
    assert (len(res)==model.nAtom)

    for n, atom in enumerate(model.atom):
        cmd_str = f"alter {sel_str} and resi {atom.resi}, b='{res[n]}'"
        if not print_cmd:
            cmd.do(cmd_str)
        else:
            print(cmd_str)
apply_erp = apply_metric_from_npz
cmd.extend("apply_erp", apply_metric_from_npz)

def load_with_error_metrics(pdb_path, metric_type='lddt', selection='all', npz_path=None, print_cmd=False):
    """Loads a PDB and the error metrics"""
    
    if npz_path is None:
        npz_path=pdb_path.replace('.pdb', '.npz')
        
    cmd_str=f"load {pdb_path}"
    if print_cmd:
        print(cmd_str)
    else:
        cmd.do(cmd_str)
    apply_metric_from_npz(npz_path, metric_type=metric_type, sel_str=selection, print_cmd=print_cmd)
        
load_erp = load_with_error_metrics
cmd.extend("load_erp", load_with_error_metrics)

r"""
dela
run Z:\projects\truncator\truncator\pymol_utils.py
load_erp N1_C3_DR64_HS57_cryoat_mALb8_cutT1_BA__c_term__o7__r34_0001_0004_0001_asym.pdb, metric_type=lddt
spectrum b, rainbow, minimum=0.2, maximum=1
print(get_alignment_map("/ZCON_1__numH3__from-22.38__to07.23__grAB-CD//A","DHR08_trim"))
"""




#taken from http://www.protein.osaka-u.ac.jp/rcsfp/supracryst/suzuki/jpxtal/Katsutani/InterfaceResidues.py
from pymol import stored
 
def interface_residues(cmpx, cA='c. A', cB='c. B', cutoff=1.0, selName="interface"):
    """
    interfaceResidues -- finds 'interface' residues between two chains in a complex.
 
    PARAMS
        cmpx
            The complex containing cA and cB
 
        cA
            The first chain in which we search for residues at an interface
            with cB
 
        cB
            The second chain in which we search for residues at an interface
            with cA
 
        cutoff
            The difference in area OVER which residues are considered
            interface residues.  Residues whose dASA from the complex to
            a single chain is greater than this cutoff are kept.  Zero
            keeps all residues.
 
        selName
            The name of the selection to return.
 
    RETURNS
        * A selection of interface residues is created and named
            depending on what you passed into selName
        * An array of values is returned where each value is:
            ( modelName, residueNumber, dASA )
 
    NOTES
        If you have two chains that are not from the same PDB that you want
        to complex together, use the create command like:
            create myComplex, pdb1WithChainA or pdb2withChainX
        then pass myComplex to this script like:
            interfaceResidues myComlpex, c. A, c. X
 
        This script calculates the area of the complex as a whole.  Then,
        it separates the two chains that you pass in through the arguments
        cA and cB, alone.  Once it has this, it calculates the difference
        and any residues ABOVE the cutoff are called interface residues.
 
    AUTHOR:
        Jason Vertrees, 2009.		
    """
    # Save user's settings, before setting dot_solvent
    oldDS = cmd.get("dot_solvent")
    cmd.set("dot_solvent", 1)
 
    # set some string names for temporary objects/selections
    tempC, selName1 = "tempComplex", selName+"1"
    chA, chB = "chA", "chB"
 
    # operate on a new object & turn off the original
    cmd.create(tempC, cmpx)
    cmd.disable(cmpx)
 
    # remove cruft and inrrelevant chains
    cmd.remove(tempC + " and not (polymer and (%s or %s))" % (cA, cB))
 
    # get the area of the complete complex
    cmd.get_area(tempC, load_b=1)
    # copy the areas from the loaded b to the q, field.
    cmd.alter(tempC, 'q=b')
 
    # extract the two chains and calc. the new area
    # note: the q fields are copied to the new objects
    # chA and chB
    cmd.extract(chA, tempC + " and (" + cA + ")")
    cmd.extract(chB, tempC + " and (" + cB + ")")
    cmd.get_area(chA, load_b=1)
    cmd.get_area(chB, load_b=1)
 
    # update the chain-only objects w/the difference
    cmd.alter( "%s or %s" % (chA,chB), "b=b-q" )
 
    # The calculations are done.  Now, all we need to
    # do is to determine which residues are over the cutoff
    # and save them.
    stored.r, rVal, seen = [], [], []
    cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')
 
    cmd.enable(cmpx)
    cmd.select(selName1, None)
    for (model,resi,diff) in stored.r:
        key=resi+"-"+model
        if abs(diff)>=float(cutoff):
            if key in seen: continue
            else: seen.append(key)
            rVal.append( (model,resi,diff) )
            # expand the selection here; I chose to iterate over stored.r instead of
            # creating one large selection b/c if there are too many residues PyMOL
            # might crash on a very large selection.  This is pretty much guaranteed
            # not to kill PyMOL; but, it might take a little longer to run.
            cmd.select( selName1, selName1 + " or (%s and i. %s)" % (model,resi))
 
    # this is how you transfer a selection to another object.
    cmd.select(selName, cmpx + " in " + selName1)
    # clean up after ourselves
    cmd.delete(selName1)
    cmd.delete(chA)
    cmd.delete(chB)
    cmd.delete(tempC)
    # show the selection
    cmd.enable(selName)
 
    # reset users settings
    cmd.set("dot_solvent", oldDS)
 
    return rVal
 
cmd.extend("interface_residues", interface_residues)


def get_chainbreaks(objname, cutoff_A=2, cmd=None):
    """Returns the resi where the chain break is. The next residue is in a different place 
    N --chain_break-- N+1  .Returns N."""
    import numpy as np
    if cmd is None:
        import pymol
        cmd = pymol.cmd
    C_atoms = cmd.get_coords(objname+" and name C", 1)
    N_atoms = cmd.get_coords(objname+" and name N", 1)
    #subastract the C from the N of the next residue
    distances = np.sum((C_atoms[:-1]-N_atoms[1:])**2, axis=1)
    #len(distances), min(distances), np.mean(distances) ,max(distances)
    breaks = distances > cutoff_A**2
    return breaks.nonzero()[0]+1


aa_3to1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

def get_seq_alignment_map(from_sele, to_sele, aln_name=None, cmd=pymol.cmd):
    """Gets the sequnce mapping from calling align. Returns a tuple of CA atoms"""
    model1_name = cmd.get_object_list(from_sele)
    assert len(model1_name)==1
    model1_name = model1_name[0]
    model2_name = cmd.get_object_list(to_sele)
    assert len(model2_name)==1   
    model2_name = model2_name[0]

    if aln_name is None:
        aln_name = f"aln__{model1_name}__{model2_name}"
    
    cmd.align(f'({to_sele}) and name CA', f'({from_sele}) and name CA', cutoff=2.5, cycles=5, object=aln_name)

    #cmd.super(f'({to_sele}) and name CA', f'({from_sele}) and name CA', object=aln_name)
    raw_aln  = cmd.get_raw_alignment(aln_name)
    res = [] 
    for idx1, idx2 in raw_aln:
        from_at = cmd.get_model(f"index {idx1[1]} and {idx1[0]}").atom[0]
        from_at.parent_model = idx1[0]
        to_at = cmd.get_model(f"index {idx2[1]} and {idx2[0]}").atom[0]
        to_at.parent_model = idx2[0]

        if   model1_name == idx1[0] and  model2_name == idx2[0]:
            res.append((from_at, to_at))
        elif model2_name == idx1[0] and  model1_name == idx2[0]:
            res.append((to_at, from_at))
        else:   
            raise "Models names don't match"
    


    return res

def print_seq_alignement_map(smap):
    """Prints a tuple of atoms as alignment map"""
    model1 = smap[0][0].parent_model 
    model2 = smap[0][1].parent_model
    
    print(f'{model1} -> {model2}')
    for at1, at2 in smap:
        print(f'{at1.chain} {at2.resi_number:3d} ({aa_3to1[at1.resn]}) -> {at2.chain} {at2.resi_number:3d} ({aa_3to1[at2.resn]})')

def get_seq_aln_map_differences(smap):
    """Returns only the tuples that have a diffrent resisue name (resn) """
    res=[]
    for at1, at2 in smap:
        if at1.resn != at2.resn:
            res.append((at1, at2))
    return res
        


def get_selections_from_seq_aln_map(smap):
    """Returns the from and to selections strings (based on atom index and name)"""
    sel1_idx = []
    sel2_idx = []
    for at1, at2 in smap:
        sel1_idx.append(str(at1.index))
        sel2_idx.append(str(at2.index))
    sel1_idx_str = '+'.join(sel1_idx)
    sel2_idx_str = '+'.join(sel2_idx)
    
    #each selection can only be from one model, so only look at first atom
    #TODO rewrite using get_object_list
    model1 = smap[0][0].parent_model 
    model2 = smap[0][1].parent_model
    #TODO: index is not really the best selection since it changes on sort and mutations 
    sel1 = f'byresi(index {sel1_idx_str} and {model1})' 
    sel2 = f'byresi(index {sel2_idx_str} and {model2})' 
    
    return sel1, sel2
    
def filter_seq_aln_map_on_selection(smap, sel, smap_index=0):
    """Only keep those alignments that are both in sel1 and sel2"""
    sel_model = cmd.get_model(f"{sel} and name CA")
    sel_indexes = set([atom.index for atom in sel_model.atom])
    res = []
    for link in smap:    
        if link[smap_index].index in sel_indexes:
            res.append(link)
    
    return res

def show_align_diff(sel1, sel2, col1='same', col2='same', make_sel=True, print_sel=True, cmd=pymol.cmd):
    """Finds the diffrences in sequence (after alignment) and colors them if color is given"""
    smap = get_seq_alignment_map(sel1, sel2) 
    #print_seq_alignement_map(smap)

    smap = get_seq_aln_map_differences(smap)
    

    if len(smap) == 0:
        print('NO DIFFERENCES')
        return False
    diff_sel1, diff_sel2 = get_selections_from_seq_aln_map(smap)
    
    model1 = smap[0][0].parent_model 
    model2 = smap[0][1].parent_model
    
    if make_sel:
        cmd.select(f'{model1}_diff_{model2}', diff_sel1)
        cmd.select(f'{model2}_diff_{model1}', diff_sel2)
        
    if not col1 is None:
        if col1 != 'same':
            cmd.color(col1, diff_sel1)
        cmd.show('licorice', diff_sel1)
        cmd.hide('licorice', f" ({diff_sel1}) and (h. and (e. c extend 1))")
    if not col2 is None:
        if col2 != 'same':
            cmd.color(col2, diff_sel2)
        cmd.show('licorice', diff_sel2)
        cmd.hide('licorice', f" ({diff_sel2}) and (h. and (e. c extend 1))")
    
    cmd.do('util.cnc')

    if print_sel:
        print_seq_alignement_map(smap)        
    return True

cmd.extend("show_align_diff", show_align_diff)

def mutate_residue(sele, target_3resname, cmd=pymol.cmd):
    """Mutates the sele to target_3resname"""
    #Open wizard
    cmd.do("wizard mutagenesis; refresh_wizard") 
    
    cmd.select('sele', f'byresi ({sele})')
    cmd.get_wizard().do_select('''sele''')
    cmd.get_wizard().do_state(8)
    cmd.get_wizard().set_mode(target_3resname)
    cmd.get_wizard().do_state(1)
    
    
    cmd.get_wizard().apply()
    #TODO: Copy chi angles here
    #close wizard
    cmd.set_wizard()

cmd.extend("mutate_residue", mutate_residue)

def pymol_display(cmd=pymol.cmd):
    """Displays the pymol session in a Jupyter notebook"""
    
    from IPython.display import Image
    from tempfile import mktemp
    import os
    image_name = mktemp()+'.png'
    cmd.png(image_name, ray=0)
    im = Image(filename=image_name)
    os.remove(image_name)
    return im