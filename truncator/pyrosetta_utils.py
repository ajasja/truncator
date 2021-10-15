# print("IMPORTING pyrosetta")
#import pyrosetta
#from pyrosetta import rosetta
#import pyrosetta.distributed.io

def get_layer_selector_form_XML(tag_str):
    """Given a layer selector XML returns the layer selector object"""
    tag = pyrosetta.rosetta.utility.tag.Tag()
    tag.read(tag_str)
    datamap = pyrosetta.rosetta.basic.datacache.DataMap()
    layer_selector = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    layer_selector.parse_my_tag(tag, datamap)
    return layer_selector


def get_score_function_from_XML(snippet):
    """Given an XML snippet <ScoreFunction... a score function is returned"""
    sfl = rosetta.protocols.parser.ScoreFunctionLoader()
    tag = rosetta.utility.tag.Tag()
    tag.read(snippet)
    return sfl.create_scorefxn_from_tag(tag)


def get_layers(layer_selector, pose):
    """Returns the core, boundry and surfacace subsets (selections) for a pose"""
    layers = dict(core=[1, 0, 0], boundary=[0, 1, 0], surface=[0, 0, 1])
    result = {}
    for layer in layers.keys():
        layer_selector.set_layers(*layers[layer])
        result[layer] = layer_selector.apply(pose)

    return result


def get_residue_layer(res_id, layers_vectors):
    """Given a rosetta resid retunrs if the residue is in the surface, core or boundry layer"""
    for layer in layers_vectors.keys():
        if layers_vectors[layer][res_id] > 0:
            return layer


def visualise_selector(
    pose,
    layer_selector,
    selection_name="rosetta_sele",
    color="red",
    base_color="gray",
    base_name="all",
    do_print=True,
):
    selected_pymol = (
        pyrosetta.rosetta.core.simple_metrics.metrics.SelectedResiduesPyMOLMetric()
    )
    selected_pymol.set_residue_selector(layer_selector)
    sel_str = selected_pymol.calculate(pose).replace("rosetta_sele", selection_name)
    col_str_base = f"color {base_color}, {base_name}"
    col_str_sele = f"color red, {selection_name}"
    
    result_str = "\n".join([sel_str, col_str_base, col_str_sele])
    
    if do_print:
        print(result_str)
    
    return(result_str)


def get_selector_resids(selector, pose):
    """Takes a residue selector and returns a comma seperated string of residues IDs"""
    res_ids = rosetta.core.select.get_residues_from_subset(selector.apply(pose))
    res_ids = [str(id_) for id_ in res_ids]
    return ",".join(res_ids)


def pose_to_pdb_file(packed_pose, filename):
    with open(filename, "w+") as opdb:
        opdb.write(pyrosetta.distributed.io.to_pdbstring(packed_pose))

def add_labels_to_pose(pose, resids, label, error_on_out_of_bounds_index=False):
    pdb_info = pose.pdb_info()
    N = len(pose.residues)
    for resid in resids:
        resid = int(resid)
        if (resid<=N) or error_on_out_of_bounds_index: 
            pdb_info.add_reslabel(resid, label)
        

def add_labels_to_pdb(pdb, resids, label, out_name=None, error_on_out_of_bounds_index=False):
    pose = pyrosetta.io.pose_from_pdb(pdb)
    add_labels_to_pose(pose, resids, label, error_on_out_of_bounds_index)
    
    if out_name is None:
        out_name = pdb
    
    pose_to_pdb_file(pose, out_name)
    return pose


def add_seq_to_labels(pdb_file, out_name=None):
    """
    Saves the sequence information into pdb_labels. Labels are prefixed using res__.
    For example if the first residue is an ALA, then a label res__ALA gets added.
    """
    if out_name is None:
        out_name = pdb_file
    pose = pyrosetta.pose_from_pdb(pdb_file)
    
    
    for i in range(1, len(pose)+1):
        res_name = pose.residues[i].name()
        pose.pdb_info().add_reslabel(i, f'res__{res_name}')

    print('Saving to: '+out_name)
    pose.dump_pdb(out_name)
    return pose

def replace_res(pose_source, pose_target, res, verbose=True, orient_backbone=True, transfer_label='transfered'):
    """Copies residue `res` from `source_pose` to `target_pose`. Also copies over the labels."""
    if verbose:
        print(f"Replacing res {res} source: {pose_source.residue(res).name1()} -> target: {pose_target.residue(res).name1()}")  
    pose_target.replace_residue(res, pose_source.residue(res), orient_backbone=orient_backbone)    
    #copy over tha labels
    source_lab = pose_source.pdb_info().get_reslabels(res)
    tpi = pose_target.pdb_info()
    tpi.clear_reslabel(res)
    for lab in source_lab:
        tpi.add_reslabel(res, lab)
    if transfer_label:
        tpi.add_reslabel(res, transfer_label)