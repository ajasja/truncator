"""Main high-level truncator module"""
import truncator

def score_interface(input_pdb, output_score_file=None, script_name='truncator/xml/03_FD_FR_and_scoring_benchmark_v3_digs_updated_schema_BUNS3_IA_quick.xml', 
                                                  rosetta_bin="/software/rosetta/latest/bin/rosetta_scripts",
                                                  dont_cleanup=False,
                                                  out_dir=None):
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
        output_score_file = os.path.splitext(input_pdb)[0]+"_IA.sc"
        
    res = truncator.run_rosetta_script(rosetta_bin, script_name, input_pdb, tmp_dir, 
                       extra_flags="-holes:dalphaball truncator/scripts/DAlphaBall.gcc")
    
    score_file=truncator.add_to_score_from_log(log_file=res.log_file, score_file=res.score_file, out_file=output_score_file)
    if not dont_cleanup:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    return score_file
    
    
