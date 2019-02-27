"""Utility functions for truncator"""
import os, errno
import pandas as pd
import numpy as np
from glob import glob
import truncator

# taken from http://stackoverflow.com/a/11301781/952600
try:
    basestring  # attempt to evaluate basestring
    def is_str(s):
        return isinstance(s, basestring)
except NameError:
    def is_str(s):
        return isinstance(s, str)



def notebook_setup(init_flags=str):
    import pyrosetta
    flags = " ".join(init_flags.split())
    pyrosetta.init(flags)


def init_env(flags, workers=10, 
            init_flags_script = '/home/klimaj/ipython_notebooks/PyRosettascripts_Demo_Pre-Rosettacon_2018/pyrosettascripts_demo/data_generation/click/worker_setup.py',
            memory="4GB",
            queue="interactive"):
    """Function initializes current jupyter notebook environment with string `flags`,
    also initializing each of integer `workers` workers on a SLURM scheduled computer
    cluster with the same `flags`. 
    """
    from dask_jobqueue import SLURMCluster
    from dask.distributed import Client, progress
    flags_str = " ".join(flags.replace('\n', ' ').split())
    
    
    #%run {init_flags_script}
    notebook_setup(flags_str)
    
    cluster = SLURMCluster(cores=1, 
                           processes=1,
                           memory=memory,
                           queue=queue,
                           extra=" --preload {} ' {}'".format(init_flags_script, flags_str))
    print(cluster.job_script())
    cluster.scale(int(workers))
    client = Client(cluster)
    return client



def read_file(filename):
    with open(filename, 'r') as myfile:
        data = myfile.read()
    return data         

def read_file_lines(filename, trim=False, skip_comments=False):
    with open(filename, 'r') as myfile:
        lines = myfile.readlines()
    if trim or skip_comments:
        lines = [line.strip() for line in lines]
    if skip_comments:
        lines = [line for line in lines if not (line=='' or line.startswith('#'))]
    return lines 

def write_file(filename, str):
    with open(filename, 'w') as myfile:
        myfile.write(str)

def write_json(filename, data):
    import json
    with open(filename, 'w') as outfile:
        json.dump(data, outfile,  indent=4)

def read_json(filename):
    import json
    with open(filename, 'r') as _file:
        return json.load(_file)

def read_info_file(file_name, ext=None):
    if not ext is None:
        file_name = file_name.replace(".gz","").replace(".pdb","")
        file_name = file_name+ext
    if os.path.exists(file_name):
        return read_json(file_name)
    else:
        print(f"WARNING: info file not found: {file_name}")
        return {}

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

import os.path
def ensure_extension(file_name, ext):
    """Appends filname with ext, if there is no extension. 
    ext can contain a dot or not"""
 
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
    try:
        os.makedirs(dir_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

from itertools import chain
def flatten_list(_list):
    """Flattens a list of lists"""
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


def get_full_file_from_name(name, dir, exts=".pdb .pdb.gz".split()):
    """Tries to get the file from the name, the directory and a list of extensions. The first match that exists is returned"""
    import os
    for ext in exts:
        fname = os.path.join(dir, name)+ext
        if os.path.exists(fname):
            return fname
    
    return None

import Bio
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

GLOBAL_PDBPARSER = None
GLOBAL_PPBUILDER = None


def seq_from_pdb(file_name):
    """Returns sequence from a PDB file (all chains concatanated). Accepts .gz pdbs"""
    global GLOBAL_PDBPARSER
    global GLOBAL_PPBUILDER
    import os

    try:
        if os.path.splitext(file_name)[-1] == '.gz':
            import gzip
            f = gzip.open(file_name, 'rt')
        else:
            f = open(file_name, 'r')
        if GLOBAL_PDBPARSER is None:
            GLOBAL_PDBPARSER = PDBParser()
            GLOBAL_PPBUILDER = PPBuilder()
        
        struct = GLOBAL_PDBPARSER.get_structure('', f)
        seqs = [str(pp.get_sequence()) for pp in GLOBAL_PPBUILDER.build_peptides(struct)]
        return "".join(seqs)
    finally:
        f.close()


def read_score_file(file_name, pdb_dir=None, verbose=False, load_seq=True, skiprows=1):
    """Loads a score file and PDB sequnces to a pandas dataframe"""
    import pandas as pd

    #df = pd.read_csv(file_name, skiprows=skiprows, header=0, sep=" ")
    #df.drop(columns=['SCORE:'], inplace=True, errors='ignore')
    df = pd.read_json(file_name, orient='records', lines=True)
    if pdb_dir is None:
        pdb_dir = os.path.dirname(file_name)
    full_names = []
    seqs = []
    #charges = []
    #pIs = []
    #MWs = []
  
    for n, name in enumerate(df.decoy.values):
        if verbose:
            print(name)
        full_name = get_full_file_from_name(name, pdb_dir)
        full_names.append(full_name)
        if load_seq:
            seqs.append(seq_from_pdb(full_name))
    df['full_name'] = full_names
    
    if load_seq:
        df['seq'] = seqs
    else:
        df['seq'] = " "
    return df

def read_score_files(file_names, pdb_dir=None, verbose=False, load_seq=True, skiprows=1, cache_file=None):
    dfs = []
    
    for file_name in file_names:
        dfs.append(read_score_file(file_name, pdb_dir=pdb_dir, verbose=verbose, load_seq=load_seq, skiprows=skiprows))
    dfs = pd.concat(dfs, sort=False)
    
    return dfs


def read_score_files_with_cache(file_names, cache_file, pdb_dir=None, verbose=False, load_seq=False, skiprows=1, force_reload=False):
    if os.path.exists(cache_file) and not force_reload:
        dfs = pd.read_csv(cache_file, index_col=False)
        return dfs
    
    dfs = read_score_files(file_names, pdb_dir=pdb_dir, verbose=verbose, load_seq=load_seq, skiprows=skiprows)
    
    dfs.to_csv(cache_file, index=False)
    return dfs


def find_input(query, inputs, find_all=False):
    """Finds the line containing the query string"""
    res=[]
    for inp in inputs:
        if query in inp:
            if find_all:
                res.append(inp)
            else:
                return inp
    
    if find_all:
        return res
    else:
        return None

def find_unfinished_logs(adir):
    logs = glob(adir+"/*.log")
    result = []
    for log in logs:
        log_lines = truncator.read_file_lines(log, trim=True)
        if truncator.find_input("no more batches to process...", log_lines)==None:
            result.append(log)
    return result

def find_no_output_pdbs(adir, ignore_unfinished=True):
    """Finds basenames that have a log file but no pdb result files"""
    pdbs = glob(adir+"/*.pdb")
    pdbs = [pdb.replace("_0001.pdb",".pdb").replace(".pdb","") for pdb in pdbs]
    #print(pdbs)
    logs = glob(adir+"/*.log")
    logs = [log.replace(".log","") for log in logs]
    
    
    #print(scs)
    #find names that have a score file, but no pdb
    result = set(logs)-set(pdbs)
    if ignore_unfinished:
        ulogs = find_unfinished_logs(adir)
        ulogs = [log.replace(".log","") for log in ulogs]
        result = result - set(ulogs)
    return sorted(result)

def clean_unfinished_logs(adir):
    """Cleans unfinished logs"""
    logs = truncator.find_unfinished_logs(adir)
    import os
    for log in logs:
        os.remove(log)

def lines_grep(pattern, lines):
    import re
    reg = re.compile(pattern)
    return [l for l in lines if reg.match(l)]

def pp_flags(cmd):
    print(cmd.replace(" -", " \\\n-"))