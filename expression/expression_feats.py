from pyChimera_main.chimera import *
import codonbias as cb
import pandas as pd
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta
import re
from time import time
import pickle
import numpy as np
from itertools import product
from collections import Counter
import logging
from datetime import datetime
import os
import warnings
warnings.filterwarnings("ignore")
from Bio import BiopythonWarning
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

folder = BASE_DIR

def find_overlap_genes(cdss,start,stop,chrm):
    start_matches = (cdss.start<=start) & (start<=cdss.stop) & (chrm==cdss.chrm)
    stop_matches = (cdss.start<=stop) & (stop<=cdss.stop) & (chrm==cdss.chrm)
    df_match = cdss[start_matches | stop_matches]
    return df_match

def get_codons(df_match,start,stop,w):
    for i in df_match.index:
        chrm,strand_orf_start,orf_stop = df_match.loc[i,['chrm','strand','start','stop']]
        overlap_start = max(start,df_match)
        # seq = chrs[chrm][]

#from codonbias module
def geomean(seq,aa_or_codon,freqs):
    """
    Compute the geometric mean based on codon scores given in
    `log_weights` (weights in logarithmic scale), and codon counts give
    in `counts`.

    input:
    seq : string (aa or codons)
    aa_or_codons : bool
    freqs : dict of freqs for aa/codons
    
    to be calc (orig parameters)
    ----------
    log_weights : pandas.Series
        Codon scores in logarithmic scale, with codons as index and scores
        as values.
    counts : pandas.Series
        Codon counts, with codons as index and counts as values.

    Returns
    -------
    float
        Geometric mean.
    """
    if aa_or_codon=='aa':
        seq = seq.replace('*', '')
        if len(seq)==0:
            return 0
        curr_idxs = list(set(seq))
        log_weights = pd.Series(np.log([freqs[x] for x in set(seq)]), index=curr_idxs)
        counts = pd.Series(Counter(seq))
    elif aa_or_codon=='codon':
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
        curr_idxs = list(set(codons))
        log_weights = pd.Series(np.log([freqs[x] for x in set(codons)]), index=curr_idxs)
        counts = pd.Series(Counter(codons))
    nn = log_weights.index[np.isfinite(log_weights)]
    return np.exp((log_weights[nn] * counts.reindex(nn)).sum() / counts.reindex(nn).sum())

def print_log(x):
    print(x)
    # logging.info(x)

def calc_chimera_CAI_scores(row,SA_cod,aa_freqs,codon_freqs,cai,genome,cdss,orfs,organism='human',calc_CAI_chimera='all'):
    # wtf_idxs=[]
    if calc_CAI_chimera=='all':
        num_feats = 24
    elif calc_CAI_chimera=='CAI':
        num_feats = 18
    elif calc_CAI_chimera=='chimera':
        num_feats = 6
    if row.name%50==0:
        print_log(row.name)
    scores = []
    g_rna_info = row['g_rna_info']
    
    if organism=='human':
        chrm_prefix = 'chr'
        split_g_rna_info = g_rna_info.split('_')
        char_to_split = '_'
        if len(split_g_rna_info)!=4 or (split_g_rna_info[0] not in [str(i) for i in range(1,24)] and split_g_rna_info[0] not in ['X','Y']):
            return [[] for _ in range(num_feats)]
    else:
        chrm_prefix = ''
        char_to_split=';'

    chrm,site_start,site_stop,site_strand = g_rna_info.split(char_to_split)
    chimera_scores = {'':0}
    codon_freqs_scores = {'':0}
    aa_freqs_scores = {'':0}
    cai_scores = {'':0}
    for w in [0,10,20]:
        # print_log('w',w)
        chimera_mean_score = []
        aa_freqs_mean_score = []
        codon_freqs_mean_score = []
        cai_mean_score = []
        cai_max_score = 0
        chimera_max_score = 0
        aa_freqs_max_score = 0
        codon_freqs_max_score = 0
        # cai_max_expr_score = 0
        # aa_freqs_expr_max_score = 0
        # codon_freqs_expr_max_score = 0
        # chimera_max_expr_score = 0
        site_start = int(site_start)-w
        site_stop = int(site_stop)+w
        df_match = find_overlap_genes(cdss,site_start,site_stop,chrm)
        for j in df_match.index:
            # if j!=1:
            #     continue
            # print_log('j',j,'w',w)
            cds_start,cds_stop,cds_strand,transcript_id = df_match.loc[j,['start', 'stop', 'strand','transcript_id']]
            curr_orf_match = orfs[orfs['id']==transcript_id]
            if curr_orf_match.shape[0]!=1:
                print('WTF')
                wtf_idxs.append([df_name,row.name,j,'curr_orf_match not 1'])
            orf_seq = curr_orf_match.seq.iloc[0]
            curr_cds_seq = str(genome[chrm_prefix+chrm][cds_start:cds_stop])
            curr_cds_coords = list(range(cds_start,cds_stop+1))
            if cds_strand=='-':
                curr_cds_seq = Seq.reverse_complement(curr_cds_seq)
                curr_cds_coords.reverse()
            overlap_start = max(site_start,cds_start)
            overlap_stop = min(site_stop,cds_stop)
            overlap_start_idx = curr_cds_coords.index(overlap_start)
            overlap_stop_idx = curr_cds_coords.index(overlap_stop)
            if overlap_start_idx>overlap_stop_idx:
                overlap_start_idx, overlap_stop_idx = overlap_stop_idx, overlap_start_idx
            overlap_seq = curr_cds_seq[overlap_start_idx:overlap_stop_idx].upper()
            start_overlap_in_orf = orf_seq.index(overlap_seq)
            start_overlap_in_orf = start_overlap_in_orf - (start_overlap_in_orf%3)
            end_overlap_in_orf = start_overlap_in_orf+len(overlap_seq)
            end_overlap_in_orf = end_overlap_in_orf -(end_overlap_in_orf%3)
            codons_seq = orf_seq[start_overlap_in_orf:end_overlap_in_orf]
            codons_seq = ''.join([codons_seq[i_c:i_c+3] for i_c in range(0,len(codons_seq),3) if 'N' not in codons_seq[i_c]])
            # print(len(codons_seq)%3==0,codons_seq)
            if len(codons_seq)==0:
                continue
            aa_seq = str(Seq.Seq(codons_seq).translate())
            if aa_seq not in str(Seq.Seq(orf_seq).translate()):
                print_log('seq not in seq')
                print_log(row.name,codons_seq,len(codons_seq)/3,start_overlap_in_orf%3)
                print_log(str(Seq.Seq(codons_seq).translate()))
                wtf_idxs.append([df_name,row.name,f'j {j}',f'w {w}''all of seq not in orf'])
            elif len(codons_seq)%3!=0:
                print_log('len not divisible by 3')
                wtf_idxs.append([df_name,row.name,f'j {j}',f'w {w}','not divisiable by 3'])
            if codons_seq in chimera_scores.keys():
                chimera_mean_score.append(chimera_scores[codons_seq])         
            else:
                # print_log('calculating...')
                # print_log(f'{df_name}, i={row.name}, j={j}, len={len(codons_seq)}, frame={start_overlap_in_orf%3}')
                # print_log(type(codons_seq),codons_seq)
                # q = str(Seq.Seq(codons_seq).translate())
                # print_log('here')
                
                target_cod = nt2codon(codons_seq)
                curr_chimera_score = calc_cARS(target_cod, SA_cod,verbose=False)
                chimera_mean_score.append(curr_chimera_score)         
                # print(curr_chimera_score)
                chimera_max_score = max(curr_chimera_score,chimera_max_score)
                chimera_scores[codons_seq] = curr_chimera_score
                curr_aa_score = geomean(aa_seq,'aa',aa_freqs)
                curr_codon_score = geomean(codons_seq,'codon',codon_freqs)
                curr_cai_score = cai.get_score(codons_seq)
                codon_freqs_mean_score.append(curr_codon_score)
                aa_freqs_mean_score.append(curr_aa_score)
                cai_mean_score.append(curr_cai_score)
                cai_max_score = max(cai_max_score,curr_cai_score)
                aa_freqs_max_score = max(curr_aa_score,aa_freqs_max_score)
                codon_freqs_max_score = max(curr_codon_score,codon_freqs_max_score)
                aa_freqs_scores[str(Seq.Seq(codons_seq).translate())] = curr_aa_score
                codon_freqs_scores[codons_seq] = curr_codon_score
                cai_scores[codons_seq] = curr_cai_score
            # break
        if cai_mean_score!=[]:
            cai_mean_score = np.mean(cai_mean_score)
        else:
            cai_mean_score=0 
            
        if aa_freqs_mean_score!=[]:
            aa_freqs_mean_score = np.mean(aa_freqs_mean_score)
        else:
            aa_freqs_mean_score=0
            
        if codon_freqs_mean_score!=[]:
            codon_freqs_mean_score = np.mean(codon_freqs_mean_score)
        else:
            codon_freqs_mean_score = 0
            
        if chimera_mean_score!=[]:
            chimera_mean_score = np.mean(chimera_mean_score)
        else:
            chimera_mean_score=0 
                        
        scores.extend([codon_freqs_mean_score,codon_freqs_max_score,
                       aa_freqs_mean_score,aa_freqs_max_score,
                       cai_mean_score, cai_max_score,chimera_mean_score,chimera_max_score])
    return scores

cols = [f'{name}_{f}{m}_{win}'
        for win in [0,10,20] 
        for name,f in zip(['codon','aa','CAI','chimera'],['freqs_','freqs_','','']) 
        for m in ['avg','max']]

def make_CUB_df():
    df = pd.DataFrame({x:[] for x in cols})
    #     df.drop(columns=[col for col in df.columns if any(x in col for x in cols)])
    # new_cols_df = pd.DataFrame({col: col_type  for col,col_type in zip(cols,isana_cols_types)}, index=df.index)
    # df = pd.concat([df, new_cols_df], axis=1)
    return df

def load_data(organism):
    print(f'loading {organism} CUB data...')
    char_to_split=';'
    if organism=='human':
        char_to_split='_'
        orfs_file_name = "Homo_sapiens.GRCh38.cds.all.fa"
        transcript_pattern = r'transcript_id\s+"(ENST\d+)"'
    elif organism=='tomato':
        orfs_file_name = "SollycM82_genes_v1.1.1.CDS.fasta"
        transcript_pattern = r'Parent=mRNA:([^;]+)'
    elif organism=='fly' or organism=='prawn':
        orfs_file_name = "cds_from_genomic.fna"    
        transcript_pattern = r'protein_id=([^;]+)'

    org_data_folder = Path(__file__).resolve().parent.parent / 'data' / f'{organism}_data'
    orfs = Fasta(org_data_folder / orfs_file_name)
    if organism=='human':
        orfs = pd.DataFrame({'seq':[str(orfs[x]) for x in list(orfs.records)], 'id':[x[:x.index('.')] for x in list(orfs.records)]})
    elif organism=='tomato':
        orfs = pd.DataFrame({'seq':[str(orfs[x]) for x in list(orfs.records)], 'id':[x[5:] for x in list(orfs.records)]})
    elif organism=='fly'  or organism=='prawn':
        orfs = pd.DataFrame({'seq':[str(orfs[x]) for x in list(orfs.records)],'id':[re.findall(r'(?:XP|YP)_\d+\.\d+', x)[0] for x in list(orfs.records)]})
    
    cai = cb.scores.CodonAdaptationIndex(ref_seq=orfs.seq.to_list())
    cdss = pd.read_csv(org_data_folder / f'{organism}_cdss.csv',dtype={'chrm': str})
    cdss['transcript_id'] = cdss['description'].str.extract(transcript_pattern)
    with open(org_data_folder / f'{organism}_codon_freqs.pkl', 'rb') as file:
        codon_freqs = pickle.load(file)
    with open(org_data_folder / f'{organism}_aa_freqs.pkl', 'rb') as file:
        aa_freqs = pickle.load(file)
    with open(org_data_folder / f'{organism}_SA_cod.pkl', 'rb') as file:
        SA_cod = pickle.load(file)
    print(f'done loading {organism} CUB related data')
    return orfs, cai, cdss, codon_freqs, aa_freqs, SA_cod, char_to_split

def calc_CUB_feats(df,genome,organism='',calc_CAI_chimera='all'):
    if organism=='':
        print('no organism chosen!','pick an organism from: human, fly, prawn, tomato')
        return df
    orfs, cai, cdss, codon_freqs, aa_freqs, SA_cod, char_to_split = load_data(organism)
    print('started calculating CUB feats!')
    output_df = df.apply(lambda row:
calc_chimera_CAI_scores(row,SA_cod,aa_freqs,codon_freqs,cai,
                        genome,cdss,orfs,
                        organism=organism,calc_CAI_chimera=calc_CAI_chimera),
                      axis=1, result_type='expand')
    output_df.columns = cols
    print('done calculating CUB feats!')

    return output_df

