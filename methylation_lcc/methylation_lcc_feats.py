from Bio import motifs
import Bio.SeqUtils.lcc as lcc
import numpy as np
import pandas as pd
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

base_pairs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def reverse_complement(x):
    return str(''.join(base_pairs.get(base, base) for base in reversed(x)))

with open(BASE_DIR / 'topEnriched.313.meme.txt', 'r') as handle:
    motifs_methylation = motifs.parse(handle, "minimal")

pwm_methylation = [x.counts.normalize(pseudocounts={'A':0.295, 'C': 0.205, 'G': 0.205, 'T': 0.295}) for x in motifs_methylation]
pssm_methylation = [x.log_odds(motifs_methylation.background) for x in pwm_methylation]


def calc_methyl_vector(motifs_methylation, row):
    seq = row.loc['target_seq'][:20]
    pre_seq = row.loc['up_seq']
    post_seq = row.loc['down_seq']

    if type(seq) == float:
        return row

    len_seq = len(seq)

    forward_seq = seq+post_seq[:15]
    backward_seq = reverse_complement(seq)+ reverse_complement(pre_seq[-15:])


    results_arr = motifs_methylation[0].pssm.calculate(backward_seq)[:len_seq].reshape(1, -1)
    for motif in motifs_methylation:
        np.append(results_arr, motif.pssm.calculate(backward_seq)[:len_seq].reshape(1, -1), axis=0)

    backward_vec = np.fliplr(np.amax(results_arr, axis=0).reshape(1, -1))


    results_arr = motifs_methylation[0].pssm.calculate(forward_seq)[:len_seq].reshape(1, -1)
    for motif in motifs_methylation:
        np.append(results_arr, motif.pssm.calculate(forward_seq)[:len_seq].reshape(1, -1), axis =0)

    np.append(results_arr, backward_vec, axis =0)

    max_vec = np.amax(results_arr, axis=0)

    for ii in range(len_seq):
        row[f'methylation_{ii}'] = round(max_vec[ii], 4)

    row[f'methylation_max'] = round(np.max(max_vec), 4)
    row[f'methylation_mean'] = round(np.mean(max_vec), 4)
    row[f'methylation_std'] = round(np.std(max_vec), 4)

    return row

def lcc_features(row):
    if type(row.loc['target_seq']) == float:
        return row
    curr_seq = row.loc['up_seq'][-10:]+row.loc['target_seq'][:20]+row.loc['down_seq'][:10]
    lcc_features = lcc.lcc_mult(curr_seq, 10)

    for ii, curr_lcc in enumerate(lcc_features):
        row[f'LCC_{ii}'] = round(curr_lcc, 4)
    return row

def calc_methylation_feats(df):
    df = df.apply(lambda row: calc_methyl_vector(motifs_methylation, row), axis = 1)
    return df
def calc_lcc_feats(df):
    df = df.apply(lambda row: lcc_features(row), axis = 1)
    return df

