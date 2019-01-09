"""Utility functions for truncator"""
import os
import pandas as pd
import numpy as np


def read_file(filename):
    with open(filename, 'r') as myfile:
        data = myfile.read()
    return data         

def read_file_lines(filename):
    with open(filename, 'r') as myfile:
        lines = myfile.readlines()
    return lines 

def write_file(filename, str):
    with open(filename, 'w') as myfile:
        myfile.write(str)

def write_json(filename, data):
    import json
    with open(filename, 'w') as outfile:
        json.dump(data, outfile,  indent=4)


def replace_extension(file_name, ext):
    """Replace extension of file_name with ext. Ext can contain a dot or not"""
    import os.path
    base = os.path.splitext(file_name)[0]

    if len(ext) == 0:
        dot = ""
    else:
        if ext[0] == ".":
            dot = ""
        else:
            dot = "."

    return base + dot + ext


def ensure_extension(file_name, ext):
    """Appends filname with ext, if there is no extension. 
    ext can contain a dot or not"""
    import os.path
    base, existing_ext = os.path.splitext(file_name)

    if existing_ext != '':
        return file_name

    if ext[0] == ".":
        dot = ""
    else:
        dot = "."

    return base + dot + ext


def ensure_file_name(file_name, base_file_name, ext):
    """If file name is None then returns base file name with replaced extesnion. 
    Otherise it calls ensure extension on the file name. Ext can contain a dot or not"""
    if file_name is None:
        return replace_extension(base_file_name, ext)
    else:
        return ensure_extension(file_name, ext)


def basename_noext(filename):
    """Returns the basename without extension"""
    base_name = os.path.basename(filename) 
    base, existing_ext = os.path.splitext(base_name)
    return base

def make_dirs(dir_path):
    """Makes missing subdirs if they don't exists."""
    import os, errno
    try:
        os.makedirs(dir_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def flatten_list(_list):
    """Flattens a list of lists"""
    from itertools import chain
    return list(chain.from_iterable(_list))


import numpy as np
def find_helices(sec_struct):
    """Returns 0 based python slice convention"""

    bits = np.array(sec_struct, dtype='S1') == b'H'
    # make sure all runs of ones are well-bounded
    bounded = np.hstack(([0], bits, [0]))
    # get 1 at run starts and -1 at run ends
    difs = np.diff(bounded)
    run_starts = np.where(difs > 0)[0]
    run_ends = np.where(difs < 0)[0] 
    return np.dstack((run_starts,run_ends))[0]

def python_to_pymol_range(_from,_to):
    _from=_from+1
    return (_from,_to)
def list_to_pymol_sel_str(_list, join_char="+"):
    _list_str = [str(ele) for ele in _list]
    return join_char.join(_list_str)
#list_to_pymol_sel_str([1,2,3])=="1+2+3"