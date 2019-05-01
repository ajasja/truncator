"""Main high-level truncator module"""
import truncator
import os
import tempfile
import shutil


def score_interface(input_pdb, output_score_file=None, script_name='truncator/xml/03_FD_FR_and_scoring_benchmark_v3_digs_updated_schema_BUNS3_IA_quick.xml', 
                                                  rosetta_bin="/software/rosetta/latest/bin/rosetta_scripts",
                                                  dont_cleanup=False,
                                                  out_dir=None,
                                                  tee=False, test_run=False, skip_existing=False):
    """Score interface metrics for PDB"""   

    if out_dir is None:
        tmp_dir = tempfile.mkdtemp()
    else:
        tmp_dir = out_dir
        dont_cleanup = True

    truncator.make_dirs(tmp_dir)


    if output_score_file is None:
        if out_dir is None:
            output_score_file = os.path.abspath(os.path.splitext(input_pdb)[0]+"_IA.sc")
        else:
            output_score_file = os.path.abspath(out_dir + "/"+ truncator.basename_noext(input_pdb)+"_IA.sc")
        
    if skip_existing and os.path.exists(output_score_file):
        print(f"Output score file exists: {output_score_file}")
        return output_score_file

    res = truncator.run_rosetta_script(rosetta_bin, script_name, input_pdb, tmp_dir, 
                       extra_flags="-holes:dalphaball truncator/scripts/DAlphaBall.gcc -out:no_nstruct_label", tee=tee, test_run=test_run)
    
    score_file=truncator.add_to_score_from_log(log_file=res.log_file, score_file=res.score_file, out_file=output_score_file)
    
    if not dont_cleanup:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    return score_file
    


def reloop(input_pdb, out_dir=None, script_name='truncator/xml/12_redesign_loops_tl_gen_profile_SymAnn.xml', 
                                                  rosetta_bin="/software/rosetta/latest/bin/rosetta_scripts",
                                                  extra_flags="-chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm",
                                                  structure_store="/home/brunette/DBs/hdf5/ss_grouped_vall_helix_shortLoop.h",
                                                  log_file=None,
                                                  score_file=None,
                                                  dont_cleanup=False,
                                                  tee=False, test_run=False, verbose=True, skip_existing=False,
                                                  add_suffix = True, clean_existing=True,
                                                  chain_connections='[A+B,C+D],[B+A,C+D],[A+B,D+C],[B+A,D+C]',
                                                  loopLengthRange="2,5",
                                                  resAdjustmentRangeSide1="-1,1", 
                                                  resAdjustmentRangeSide2="-1,1",
                                                  allowed_loop_abegos="AGBA,ABBA,AGBBA,ABABA,ABBBA,AGABBA,ABBBBA,AGBBBA",
                                                  RMSthreshold="0.5"):
    """Reloop files"""   

    rosetta_bin=os.path.abspath(rosetta_bin)
    script_name=os.path.abspath(script_name)
    structure_store=os.path.abspath(structure_store)
    input_pdb=os.path.abspath(input_pdb)

    cc_name = truncator.chain_connections_to_safe_name(chain_connections)
    if add_suffix:
        suffix = f"__cc{cc_name}__llr{loopLengthRange}__ar1{resAdjustmentRangeSide1}__ar2{resAdjustmentRangeSide2}__rms{RMSthreshold}"
    else:
        suffix = ""


    base_name = truncator.basename_noext(input_pdb) 
    out_dir = out_dir+"/"+base_name+suffix
    truncator.make_dirs(out_dir)
    
    if score_file is None: score_file =  base_name+suffix+".sc"
    if log_file is None: log_file =      base_name+suffix+".log"
    pdb_file   =                         base_name+suffix+".pdb"

    

    if clean_existing:
        truncator.remove_file(out_dir+"/"+score_file)
        truncator.remove_file(out_dir+"/"+log_file)
        truncator.remove_file(out_dir+"/"+pdb_file)

    if tee: 
        redir_str = '| tee'
    else:
        redir_str = '>'

    script_vars="chain_connections loopLengthRange resAdjustmentRangeSide1 resAdjustmentRangeSide2 allowed_loop_abegos RMSthreshold".split()

    script_vars_str = ''
    vals = locals()
    for sv in script_vars:
        script_vars_str += f"{sv}=\"{vals[sv]}\" "


#   -out:path:pdb {out_dir} \

    cmd = f" \
    -parser:protocol {script_name} -s {input_pdb} \
    -indexed_structure_store:fragment_store  {structure_store} \
    -out:suffix {suffix} \
    -out:file:scorefile  {score_file} \
    -beta \
    -in:file:fullatom \
    -renumber_pdb 1 \
    -out:file:pdb_comments true \
    -run:preserve_header true \
    -out:file:scorefile_format json \
    -out:pdb \
    -parser:script_vars {script_vars_str} \
    -out:no_nstruct_label \
    {extra_flags} \
    {redir_str} {log_file}".replace('    ','')
    

        

    cmd = f"cd {out_dir}; {rosetta_bin} {cmd}"


    if verbose:
       truncator.pp_flags(cmd)

    if skip_existing and os.path.exists(log_file):
        print("SKIPPING: "+log_file)
        
    if not test_run:
        status = os.system(cmd)
    else:
        status = 0
    

    return truncator.ScriptRunResult(log_file, score_file, None, status, cmd)
    #return log_file    


def fix_surface(input_pdb, out_dir=None, script_name='truncator/xml/31_fix_surf.xml', 
                                                  rosetta_bin="/software/rosetta/latest/bin/rosetta_scripts",
                                                  extra_flags="-chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm",
                                                  log_file=None,
                                                  score_file=None,
                                                  out_pdb=None,
                                                  tee=False, test_run=False, verbose=True, skip_existing=False, clean_existing=True,
                                                  add_suffix = True,
                                                  charge_constraints_chA = "",
                                                  charge_constraints_chB = "",
):
    """Fix the surface output"""   

    rosetta_bin=os.path.abspath(rosetta_bin)
    script_name=os.path.abspath(script_name)
    input_pdb=os.path.abspath(input_pdb)
    charge_constraints_chA=os.path.abspath(charge_constraints_chA)
    charge_constraints_chB=os.path.abspath(charge_constraints_chB)

    if add_suffix:
        suffix = f"__ch-5"
    else:
        suffix = ""


    base_name = truncator.basename_noext(input_pdb) 
    out_dir = out_dir+"/"+base_name+suffix
    truncator.make_dirs(out_dir)
    
    if score_file is None: score_file =  base_name+suffix+".sc"
    if log_file is None: log_file =      base_name+suffix+".log"
    pdb_file   =                         base_name+suffix+".pdb"
  

    if clean_existing:
        truncator.remove_file(out_dir+"/"+score_file)
        truncator.remove_file(out_dir+"/"+log_file)
        truncator.remove_file(out_dir+"/"+pdb_file)

    if tee: 
        redir_str = '| tee'
    else:
        redir_str = '>'

    script_vars="charge_constraints_chA charge_constraints_chB".split()

    script_vars_str = ''
    vals = locals()
    for sv in script_vars:
        script_vars_str += f"{sv}=\"{vals[sv]}\" "
#    -parser:script_vars {script_vars_str} \

#   -out:path:pdb {out_dir} \

    cmd = f" \
    -parser:protocol {script_name} -s {input_pdb} \
    -out:suffix {suffix} \
    -out:file:scorefile  {score_file} \
    -beta \
    -in:file:fullatom \
    -renumber_pdb 1 \
    -out:file:pdb_comments true \
    -run:preserve_header true \
    -out:file:scorefile_format json \
    -out:pdb \
    -out:no_nstruct_label \
    -parser:script_vars {script_vars_str} \
    {extra_flags} \
    {redir_str} {log_file}".replace('    ','')
    

    cmd = f"cd {out_dir}; {rosetta_bin} {cmd}"


    if verbose:
       truncator.pp_flags(cmd)

    if skip_existing and os.path.exists(log_file):
        print("SKIPPING: "+log_file)
        
    if not test_run:
        status = os.system(cmd)
    else:
        status = 0
    

    return truncator.ScriptRunResult(log_file, score_file, None, status, cmd)
    #return log_file    