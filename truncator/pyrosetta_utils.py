# print("IMPORTING pyrosetta")
import pyrosetta
from pyrosetta import rosetta


def get_layer_selector_form_XML(tag_str):
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
    layers = dict(core=[1, 0, 0], boundary=[0, 1, 0], surface=[0, 0, 1])
    result = {}
    for layer in layers.keys():
        layer_selector.set_layers(*layers[layer])
        result[layer] = layer_selector.apply(pose)

    return result


def get_residue_layer(res_id, layers_vectors):
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
):
    selected_pymol = (
        pyrosetta.rosetta.core.simple_metrics.metrics.SelectedResiduesPyMOLMetric()
    )
    selected_pymol.set_residue_selector(layer_selector)
    print(selected_pymol.calculate(pose).replace("rosetta_sele", selection_name))
    print(f"color {base_color}, {base_name}")
    print(f"color red, {selection_name}")


def get_selector_resids(selector, pose):
    """Takes a residue selector and returns a comma seperated string of residues IDs"""
    res_ids = rosetta.core.select.get_residues_from_subset(selector.apply(pose))
    res_ids = [str(id_) for id_ in res_ids]
    return ",".join(res_ids)

