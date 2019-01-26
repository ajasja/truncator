"""Main high-level truncator module"""
import truncator
import os

def score_interface(input_pdb, output_score_file=None, script_name='truncator/xml/03_FD_FR_and_scoring_benchmark_v3_digs_updated_schema_BUNS3_IA_quick.xml', 
                                                  rosetta_bin="/software/rosetta/latest/bin/rosetta_scripts",
                                                  dont_cleanup=False,
                                                  out_dir=None,
                                                  tee=False, test_run=False):
    """Score interface metrics for PDB"""   
    import tempfile
    import shutil
    import os
    if out_dir is None:
        tmp_dir = tempfile.mkdtemp()
    else:
        tmp_dir = out_dir
        dont_cleanup = True
        
    if output_score_file is None:
        if out_dir is None:
            output_score_file = os.path.splitext(input_pdb)[0]+"_IA.sc"
        else:
            output_score_file = out_dir + "/"+ truncator.basename_noext(input_pdb)+"_IA.sc"
        
    res = truncator.run_rosetta_script(rosetta_bin, script_name, input_pdb, tmp_dir, 
                       extra_flags="-holes:dalphaball truncator/scripts/DAlphaBall.gcc -out:no_nstruct_label", tee=tee, test_run=test_run)
    
    score_file=truncator.add_to_score_from_log(log_file=res.log_file, score_file=res.score_file, out_file=output_score_file)
    if not dont_cleanup:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    return score_file
    
    
def reloop(input_pdb, out_dir=None, script_name='truncator/xml/10_redesign_loops_tl_gen_profile.xml', 
                                                  rosetta_bin="/software/rosetta/latest/bin/rosetta_scripts",
                                                  extra_flags="",
                                                  structure_store="/home/brunette/DBs/hdf5/ss_grouped_vall_helix_shortLoop.h",
                                                  log_file=None,
                                                  score_file=None,
                                                  dont_cleanup=False,
                                                  tee=False, test_run=False, verbose=True):
    """Score interface metrics for PDB"""   


    truncator.make_dirs(out_dir)
    base_name = truncator.basename_noext(input_pdb) 
    if score_file is None: score_file =  out_dir+"/"+base_name+".sc"
    if log_file is None: log_file =    out_dir+"/"+base_name+".log"
    
    if tee: 
        redir_str = '| tee'
    else:
        redir_str = '>'

    #
    cmd = f" \
    -parser:protocol {script_name} -s {input_pdb} \
    -out:path:all {out_dir} \
    -indexed_structure_store:fragment_store  {structure_store} \
    -out:file:scorefile {base_name}.sc \
    -beta \
    -in:file:fullatom \
    -renumber_pdb 1 \
    -overwrite \
    -out:file:pdb_comments \
    -run:preserve_header \
    -out:file:scorefile_format json \
    -out:pdb true \
    {extra_flags} \
    {redir_str} {log_file}"
    

    cmd = rosetta_bin + " " +cmd


    if verbose:
        print (cmd)
    if not test_run:
        status = os.system(cmd)
    else:
        status = 0
    

    return truncator.ScriptRunResult(log_file, score_file, None, status, cmd)    