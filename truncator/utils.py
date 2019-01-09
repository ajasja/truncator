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