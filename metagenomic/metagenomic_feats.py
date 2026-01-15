import pandas as pd
from Bio import Seq
import numpy as np
import collections
import csv
import itertools
import json
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

SPACERS = [
    ("all_spacers_cd95_cov100", "all_spacers"),
    ("Cas9_spacers_cd95_cov100", "Cas9_spacers"),
    ("Cas9a_spacers_cd95_cov100", "Cas9a_spacers"),
    ("Cas9a_spy_cluster_cd90_cov80_spacers_cd95_cov100",
     "Cas9a_spy_like_cluster_spacers"),
    ("Cas9a_spy_cluster_spacers_cd90_cov80", "Cas9a_spy_cluster_spacers"),
]

class Quantiler:
    RANKING = {
        "A": 1,
        "C": 2,
        "T": 4,
        "G": 3
    }
    OPPOSITES = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C"
    }

    def __init__(self, k, spacers_name, orig="", new="",CT_GA = ''):
        if CT_GA!='':
            CT_GA = 'tr_CT_GA-'
        if k==5:
            k_name = 'k_5-'
            folder = BASE_DIR / "5_mers"
        elif k==3:
            k_name = 'k_3-'
            folder = BASE_DIR / "3_mers"
        if spacers_name=='all_spacers':
            spacers_full_name = 'all_spacers_cd95_cov100'
            
        elif spacers_name=='Cas9_spacers':
            spacers_full_name = 'Cas9_spacers_cd95_cov100'
            
        elif spacers_name=='Cas9a_spacers':
                spacers_full_name = 'Cas9a_spacers_cd95_cov100'
            
        elif spacers_name=='Cas9a_spy_like_cluster_spacers':
                spacers_full_name = 'Cas9a_spy_cluster_cd90_cov80_spacers_cd95_cov100'
            
        elif spacers_name=='Cas9a_spy_cluster_spacers':
                spacers_full_name = 'Cas9a_spy_cluster_spacers_cd90_cov80'
            
        sorted_json_name = folder / f"{k_name}{CT_GA}{spacers_full_name}.quantiles.json"

        self.from_json(sorted_json_name)
        self.is_tr = bool(CT_GA)

        self.translate = str.maketrans(orig, new)
        self.k = k
        self.number_of_complements = self.get_number_of_complements()
        self.min_energies = None

    def to_json(self):
        with open(self.path.quantiles(self.k), "w") as file:
            json.dump(self.kmers, file)

    def from_json(self,json_path):
        self.kmers = {}
        with open(json_path, "r") as file:
            kmers = json.load(file)
            for kmer, info in kmers.items():
                self.kmers[kmer] = info[0], float(info[1])

    def get_number_of_complements(self):
        bases = [base for base in self.OPPOSITES.keys()]
        possibilities = ["".join(p) for p in itertools.product(bases, repeat=self.k)]
        result = set()
        for possibility in possibilities:
            result.add(self.get_correct_direction_complement(possibility.translate(self.translate), is_tr=self.is_tr))

        if self.is_tr:
            for i in result.copy():
                if i[-1::-1] in result and i[-1::-1] != i:
                    result.remove(i)

        return len(result)

    @classmethod
    def get_correct_direction_complement(cls, sequence, is_tr=False):
        length = len(sequence)
        if is_tr:
            for i in set(range(int(length / 2) + 1)):
                if cls.RANKING[sequence[i]] < cls.RANKING[sequence[length - i - 1]]:
                    return sequence

                if cls.RANKING[sequence[i]] > cls.RANKING[sequence[length - i - 1]]:
                    return sequence[-1::-1]

            return sequence

        for i in set(range(int(length / 2) + 1)):
            first_rank = cls.RANKING[sequence[i]]
            second_rank = cls.RANKING[cls.OPPOSITES[sequence[length - i - 1]]]
            if first_rank < second_rank:
                return sequence

            if first_rank > second_rank:
                return "".join([cls.OPPOSITES[sequence[length - i - 1]] for i in range(length)])

        return sequence

    def get_number_of_kmers(self, sequence):
        return len(sequence) - self.k + 1

    def count_reverse_complement_kmers(self, sequence):
        mers = collections.defaultdict(int)
        for i in range(self.get_number_of_kmers(sequence)):
            mers[self.get_correct_direction_complement(sequence[i: i + self.k], is_tr=self.is_tr)] += 1

        return mers

    def get_kmers_log_distribution(self, sequence):
        sequence_kmers = self.count_reverse_complement_kmers(sequence.translate(self.translate))
        number_of_mers_in_sequence = self.get_number_of_kmers(sequence)
        occurrences = {j: 0 for j in range(4)}
        frequencies = 0
        for kmer, count in sequence_kmers.items():
            if kmer in self.kmers:
                actual_kmer = kmer
            else:
                actual_kmer = str(Seq.Seq(kmer).reverse_complement())
            if actual_kmer in self.kmers:
                occurrences[self.kmers[kmer][0]] += count
                frequencies += self.kmers[kmer][1] * count
            else:
                print(f'WTF {kmer}')

        return [count / number_of_mers_in_sequence for quantile, count in occurrences.items()] + [round(frequencies,1)]

    def orig_get_kmers_log_distribution(self, sequence):
        sequence_kmers = self.count_reverse_complement_kmers(sequence.translate(self.translate))
        number_of_mers_in_sequence = self.get_number_of_kmers(sequence)
        occurrences = {j: 0 for j in range(4)}
        frequencies = 0
        for kmer, count in sequence_kmers.items():
            if kmer in self.kmers:
                occurrences[self.kmers[kmer][0]] += count
                frequencies += self.kmers[kmer][1] * count

            else:
                occurrences[4 - 1] += count
    
        return [count / number_of_mers_in_sequence for quantile, count in occurrences.items()] + [frequencies]#, sequence_kmers 

def calc_metagenomic_feats(df):
    print('started creating metagenomic features!')
    for kkk in  [5,3]:
        for spacer_i in range(5):
            for do_CT_GA in ['','do_CT_GA']:
                spacers = SPACERS[spacer_i][1]
                ugly_spacers = SPACERS[spacer_i][0]
                if do_CT_GA:
                    ugly_GCvAT = '-GCvsAT'
                    ugly_tr_CT_GA = '-tr_CT_GA'
                else:
                    ugly_GCvAT = ''
                    ugly_tr_CT_GA = ''
                a = Quantiler(kkk, spacers, orig="", new="",CT_GA=do_CT_GA)
                cols_to_check = [f'k{kkk}-{spacers}-Q{i}{ugly_GCvAT}' for i in range(1,5)]
                cols_to_check.append(f'k_{kkk}{ugly_tr_CT_GA}-{ugly_spacers}_log_freq_obs{ugly_GCvAT}')
                if do_CT_GA:
                    seqs = df['target_seq'].apply(lambda x: x.replace('C','G').replace('T','A')).to_list()
                else:
                    seqs = df['target_seq'].to_list()
                raw_results = [a.orig_get_kmers_log_distribution(z[:20]) for z in seqs]
                results = [[],[],[],[],[]]
                for i in range(5):
                    results[i] = [x[i] for x in raw_results]
                calc_burstein_df = pd.DataFrame({x:y for x,y in zip(cols_to_check,results)})
                df = pd.concat([df,calc_burstein_df],axis=1)
                # diff = df[cols_to_check]-orig_df[cols_to_check]
                # q=diff<0.00000001#00000001
                # if  ~np.all(q):# and df_name in ['T','','','']:
                #     print(f'diff {df_name} {cols_to_check}')
    print('done creating metagenomic features!')
    return df