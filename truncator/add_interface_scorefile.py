"""
This script reads the log file from the InterfaceAnalyzer mover xml and add the IA properties to the score file

"""

import sys
from glob import glob
import pandas as pd
import os


def get_interface_info(file_name, raise_exception=False):
    ### read in out file and find interface residues
    with open(file_name,'r') as file:
        d = {}
        for line in file:
            line=line.strip()
            #print lin:
            if 'TOTAL SASA:' in line:
                d['total_sasa_IA'] = float(line.split(" ")[-1])
                d['decoy'] = line.split(':')[0]
            if 'NUMBER OF RESIDUES:' in line:
                d['total_res_IA'] = float(line.split(" ")[-1])
            if 'AVG RESIDUE ENERGY:' in line:
                d['average_res_energy_IA'] = float(line.split(" ")[-1])
            if 'INTERFACE DELTA SASA:' in line:
                d['interface_delta_sasa_IA'] = float(line.split(" ")[-1])
            if 'INTERFACE HYDROPHOBIC SASA:' in line:
                d['interface_hydrophobic_sasa_IA'] = float(line.split(" ")[-1])
            if 'INTERFACE POLAR SASA:' in line:
                d['interface_polar_sasa_IA'] = float(line.split(" ")[-1])
            if 'CROSS-INTERFACE ENERGY' in line:
                d['cross_interface_energy_IA'] = float(line.split(" ")[-1])
            if 'SEPARATED INTERFACE ENERGY DIFFERENCE:' in line:
                d['separated_interface_energy_IA'] = float(line.split(" ")[-1])
            if 'DELTA UNSTAT HBONDS:' in line:
                d['delta_unsat_hbonds_IA'] = float(line.split(" ")[-1])
            if 'CROSS INTERFACE HBONDS:' in line:
                d['cross_interface_hbonds_IA'] = float(line.split(" ")[-1])
            if 'HBOND ENERGY:' in line:
                d['hbond_energy_IA'] = float(line.split(" ")[-1])
            if line.startswith('select'):
                d['pymol_interface_selection_IA']=line
    if d=={} and raise_exception:
        raise Exception(f"Log file {file_name} has no interface info!")
    return d


def add_to_score_from_log(dirname=None, log_file=None, score_file=None, out_file=None, add_nstruc_label=False):
    """Given a dir name the function loads a score file and a log file, takes the info from the log file and appends it to the score file
    Assumes score.sc, log.txt and score_IA.sc"""
    if log_file is None:
        log_file = dirname+"/log.txt"
    if score_file is None:
        score_file = dirname+"/score.sc"
    if out_file is None:
        base = os.path.splitext(score_file)[0]
        out_file = base+"_IA.sc"
    
    d=get_interface_info(log_file)

    if add_nstruc_label:
        d['decoy']=d['decoy']+"_0001"
    score = pd.read_json(score_file, lines=True)
    if not d=={}:
        interface = pd.DataFrame.from_records([d], index='decoy')
        #print(score_file)

        #score = pd.read_csv(score_file, low_memory=False, header=0, skiprows=1, sep = r'\s+', skipinitialspace=True, warn_bad_lines=False)
        #score.drop(columns=['SCORE:'],inplace=True, errors='ignore')

        
        ##take the last row and set index on description
        score = score.iloc[[-1]].set_index('decoy', drop=False)
        joined = score.join(interface, on="decoy")
    else:
        joined = score
    print("Saving to: "+out_file)
    joined.to_json(out_file, orient='records', lines=True)
    return out_file



if __name__ == "__main__":


    import argparse
    parser = argparse.ArgumentParser(description='Add info from log to score file')
    parser.add_argument('--dir', 							
                        help='The directory location of score.sc and log.txt')
    args = parser.parse_args()
    add_to_score_from_log(args.dir)