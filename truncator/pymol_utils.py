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
    cmd.iterate( f'{sel_str} and (name CA)', 'stored.ss_array.append(ss)')
    return "".join(stored.ss_array)

def get_selection_property(sel_str, property="resi"):
    """returns the secondary structure of selection as a string"""
    stored.seli_array = []
    cmd.iterate( f'{sel_str} and (name CA)', f'stored.seli_array.append({property})')
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
    print(onlyh)
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
    print(onlyh)
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
    if int(n)>0:
        return sel_first_helix(int(n), sele=sele, onlyh=onlyh, sele_name=sele_name)
    else:
        return sel_last_helix(-int(n), sele=sele, onlyh=onlyh, sele_name=sele_name)


cmd.extend("sel_first_helix", sel_first_helix)
cmd.extend("sel_last_helix", sel_last_helix)
cmd.extend("sel_terminal_helix", sel_terminal_helix)