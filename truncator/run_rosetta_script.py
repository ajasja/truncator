"""Python module to run rosetta scripts."""
"""May be completly obsolete"""
import truncator
import os
from collections import namedtuple
ScriptRunResult = namedtuple("ScriptRunResult", "log_file score_file pdb_file status run_command")


def run_rosetta_script(rosetta_bin, script_name, struct_name, out_dir, 
                       score_file=None, log_file=None, pdb_file=None, 
                       extra_flags=None,
                       verbose=True, test_run=False, tee=False, specify_output=True):
    
    base_name = truncator.basename_noext(struct_name) 

    truncator.make_dirs(out_dir)
    if score_file is None: score_file =  out_dir+"/"+base_name+".sc"
    if log_file is None: log_file =    out_dir+"/"+base_name+".log"
    if pdb_file is None: pdb_file =    out_dir+"/"+base_name+".pdb"
    
    if tee: 
        redir_str = '| tee'
    else:
        redir_str = '>'
    specify_output_str =f"\
    -out:file:o  {pdb_file} \
    -out:file:scorefile  {score_file} "
    cmd = f" \
    -parser:protocol {script_name} -s {struct_name} \
    -print_pymol_selection 0 \
    -beta \
    -in:file:fullatom \
    -renumber_pdb 1 \
    -overwrite \
    -out:file:scorefile_format json \
    {extra_flags} \
    {redir_str} {log_file}"
    
    if specify_output:
       cmd = specify_output_str + " " + cmd
    cmd = rosetta_bin + " " +cmd
    
    if verbose:
        print (cmd)
    if not test_run:
        status = os.system(cmd)
    else:
        status = 0
        
    return ScriptRunResult(log_file, score_file, pdb_file, status, cmd)    
