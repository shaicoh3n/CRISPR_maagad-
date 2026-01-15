import pandas as pd
import numpy as np
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta
import pickle
import re
import random
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from time import time
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

#DNA shape + enthalpy
with open(BASE_DIR / "dnaShape.pkl", "rb") as file:
    DNASHAPE_DICT = pickle.load(file)
DNA_PAIRS_THERMODYNAMICS = {"AA": 9.1, "AT": 8.6, "TA": 6.0, "CA": 5.8, "GT": 6.5, "CT": 7.8, "GA": 5.6, "CG": 11.9,
							"GC": 11.1, "GG": 11.0, "TT": 9.1, "TG": 5.8, "AC": 6.5, "AG": 7.8, "TC": 5.6, "CC": 11.0} #Breslauer et al.

def get_avg(l):
	return sum(l)/float(len(l))
def get_DNAshape_features(dna_seq):
    # """
    # :param dna_seq: sequence of nucleotides
    # :return: a dictionary with scores of rigidity for Major Groove Width (MGW), ProT (Propeller-Twist), Roll, and HelT (Helical-Twist).
    # The values are the scores for each pentamer/hexamer as computed by DNAshape (Zhou et al., doi:10.1093/nar/gkt437)
    # across the DNA sequence
    # """
    
    mgw = [None]
    roll = [None]
    prot = [None]
    helt = [None]

    for i in range(2, len(dna_seq) - 2):
        current_heptamer = dna_seq[i - 2: i + 3]
        current_heptamer = re.sub("N", random.choice(["A", "C", "G", "T"]), current_heptamer)
        current_nucleotide = DNASHAPE_DICT[current_heptamer]
        mgw += current_nucleotide["MGW"]
        roll += current_nucleotide["Roll"]
        prot += current_nucleotide["ProT"]
        helt += current_nucleotide["HelT"]

    ####ORIG#####
    # helt_modified = [helt[1]]
    helt_modified = []
    for i in range(2, len(helt), 2):
        helt_modified.append(get_avg(helt[i:i + 2]))
    # roll_modified = [roll[1]]
    roll_modified = []
    for i in range(2, len(roll), 2):
        roll_modified.append(get_avg(roll[i:i + 2]))
    return {"MGW": mgw[1:], "ProT": prot[1:], "Roll": roll_modified, "HelT": helt_modified}

#define columns:
def make_DNAshape_df(df):
    cols_to_drop = ['DNAshape','enthalpy']
    df = df.drop(columns=[col for col in df.columns if any(x in col for x in cols_to_drop)])
    for x in ['MGW', 'ProT', 'Roll', 'HelT']:
        columns_to_add = [f'DNAshape_{x}_{i}' for i in range(1,20)]+\
        [f'DNAshape_{x}_mean'] + [f'DNAshape_{x}_mean_ext'] + \
        [f'DNAshape_{x}_mean_up_ext'] + [f'DNAshape_{x}_mean_down_ext']
        new_cols_df = pd.DataFrame({col: np.nan for col in columns_to_add}, index=df.index)
        df = pd.concat([df, new_cols_df], axis=1)

    columns_to_add = [f'enthalpy_{i}' for i in range(1,23)]+\
    [f'enthalpy_mean'] + [f'enthalpy_mean_ext'] + \
    [f'enthalpy_mean_up_ext'] + [f'enthalpy_mean_down_ext']
    new_cols_df = pd.DataFrame({col: np.nan for col in columns_to_add}, index=df.index)
    df = pd.concat([df, new_cols_df], axis=1)
    return df

def gen_seq(chrm,start,end,strand,seq_type,genome,chrm_prefix=''):
    # if human might need to change chrm_prefix to 'chr'   
    if seq_type=='target':
        if strand=='+':
            start_delta = 0
            end_delta = 3
        else:
            start_delta = 3
            end_delta = 0
    elif seq_type=='extended':
        if strand=='+':
            start_delta = 100
            end_delta = 103
        else:
            start_delta = 103
            end_delta = 100
    elif seq_type=='down':
        if strand=='+':
            start_delta = 0
            end_delta = 103
        else:
            start_delta = 103
            end_delta = 0        
    elif seq_type=='up':
        if strand=='+':
            start_delta = 100
            end_delta = 0
        else:
            start_delta = 0
            end_delta = 100        
    seq = str(genome[chrm_prefix+chrm][start-start_delta:end+end_delta])
    if strand =='-':
        x = Seq.Seq(seq)
        x = x.reverse_complement()
        seq = str(x)
    seq = seq.upper()
    return seq

def no_N_seq(seq):
    num_N = sum([1 for x in seq if x=='N'])
    if num_N>3:
        seq = ''
    else:
        seq = re.sub("N", random.choice(["A", "C", "G", "T"]), seq)
    return seq

def gen_DNAshape_feats_per_site(chrm,start,end,strand, genome,chrm_prefix='', method='all'):
    site_feats = make_DNAshape_df(pd.DataFrame({}))
    target_seq = no_N_seq(gen_seq(chrm,start,end,strand,'target',genome,chrm_prefix))
    extended_seq = no_N_seq(gen_seq(chrm,start,end,strand,'extended',genome,chrm_prefix))
    up_seq = no_N_seq(gen_seq(chrm,start,end,strand,'up',genome,chrm_prefix))
    down_seq = no_N_seq(gen_seq(chrm,start,end,strand,'down',genome,chrm_prefix))
    if method=='dnashape' or method=='all':
        if target_seq!='':
            target_DNAshape = get_DNAshape_features(target_seq)
        if extended_seq!='':            
            extended_DNAshape = get_DNAshape_features(extended_seq)
        if up_seq!='':
            up_DNAshape = get_DNAshape_features(up_seq)
        if down_seq!='':
            down_DNAshape = get_DNAshape_features(down_seq)
        for dna_attr in ['MGW', 'ProT', 'Roll', 'HelT']:
            if target_seq!='':
                site_feats.loc[0,f'DNAshape_{dna_attr}_mean'] = np.mean(target_DNAshape[dna_attr])
                site_feats.loc[0,[f'DNAshape_{dna_attr}_{i}' for i in range(1,20)]] = target_DNAshape[dna_attr]
            if extended_seq!='':            
                site_feats.loc[0,f'DNAshape_{dna_attr}_mean_ext'] = np.mean(extended_DNAshape[dna_attr])
            if up_seq!='':
                site_feats.loc[0,f'DNAshape_{dna_attr}_mean_up_ext'] = np.mean(up_DNAshape[dna_attr])
            if down_seq!='':
                site_feats.loc[0,f'DNAshape_{dna_attr}_mean_down_ext'] = np.mean(down_DNAshape[dna_attr])
    if method=='enthalpy' or method=='all':
        if target_seq!='': 
            site_feats.loc[0,[f'enthalpy_{i}' for i in range(1,23)]] = [DNA_PAIRS_THERMODYNAMICS[target_seq[i-1:i+1]] for i in range(1, len(target_seq))]
            site_feats[f'enthalpy_mean'] = np.mean([DNA_PAIRS_THERMODYNAMICS[target_seq[i-1:i+1]] for i in range(1, len(target_seq))])
        if extended_seq!='':
            site_feats[f'enthalpy_mean_ext'] = np.mean([DNA_PAIRS_THERMODYNAMICS[extended_seq[i-1:i+1]] for i in range(1, len(extended_seq))])
        if up_seq!='':
            site_feats[f'enthalpy_mean_up_ext'] = np.mean([DNA_PAIRS_THERMODYNAMICS[up_seq[i-1:i+1]] for i in range(1, len(up_seq))])
        if down_seq!='':
            site_feats[f'enthalpy_mean_down_ext'] = np.mean([DNA_PAIRS_THERMODYNAMICS[down_seq[i-1:i+1]] for i in range(1, len(down_seq))])
    return site_feats


def calc_DNAshape(df,genome):
    chrm_prefix='' # if human - might need to change to 'chr'
    print('started creating DNAshape features!')        
    df = make_DNAshape_df(df)
    for i in df.index:
        chrm,start,end,strand = df.loc[i,'g_rna_info'].split(';')
        site_df = gen_DNAshape_feats_per_site(chrm,int(start),int(end),strand,genome,chrm_prefix)
        df.loc[i,site_df.columns] = site_df.loc[0,site_df.columns]
    print('done creating DNAshape features!')
    return df