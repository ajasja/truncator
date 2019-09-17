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


r"""
run Z:\projects\truncator\truncator\pymol_utils.py
print(get_alignment_map("/ZCON_1__numH3__from-22.38__to07.23__grAB-CD//A","DHR08_trim"))
"""
