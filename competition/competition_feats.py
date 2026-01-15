import functools
import re
# import edlib
# import multiprocess as mp
from pyfaidx import  Fasta
import pandas as pd
import numpy as np

from Bio import SeqIO, Seq
pd.options.display.max_rows = 999
pd.options.display.max_columns = 999
To evaluate the competition features, you need, in addition to the target sites, the genome and genes of the organism. The code below expects the following objects to exist:

# - "sites" DataFrame of the target site with columns sequence (the target site without PAM) and chr (its chromosome id that matches the chromosome column of the genome dataframe) and site_start (the start coordinate of the target site in the genome on the respective chromosome) and is_NGG_PAM, a boolean that indicates whether the target site has a PAM that follows the NGG pattern.
# - "genome" DataFrame is the host genome with the columns seq, chromosome
# - "genes" DataFrame with the columns seq and full_seq, where seq is the coding sequence and full_seq is the coding sequence including introns. We filter only one (the longest) isoform for each gene to not double count sites that appear in multiple isoforms.

# It will then add the competitor site counts to the sites dataframe.
# Below I give an example of the input data and compute the features for a random example site.

# #prawn
# prawn_folder = './data/prawn_data/'
# prawn_genome_path = prawn_folder + 'GCF_040412425.1_ASM4041242v1_genomic.fna'
# prawn_cdss_path =  prawn_folder + 'cds_from_genomic.fna'

# sites['sequence'] = sites.target_seq.apply(lambda x: x[:20])
# sites["chr"] = sites.g_rna_info.apply(lambda x: x.split(';')[0])
# sites["site_start"] = sites.g_rna_info.apply(lambda x: int(x.split(';')[1]))
# #make sure target_seq is 23 nts
# sites["is_NGG_PAM"] = sites.target_seq.apply(lambda x: bool(x[-2:] == 'GG'))
# def fasta_to_dataframe(filename):
#     seq_records = list(SeqIO.parse(filename, "fasta"))
#     return pd.DataFrame([{"seq": str(s.seq), "id": s.id, "name": s.name, "description": s.description} for s in seq_records])
#     genome = fasta_to_dataframe(prawn_genome_path)
#     genome["seq"] = genome["seq"].str.upper()
# q=genome['description'].str.extract(r'^(.*?M)')#
# genome["chromosome"] = q[0].apply(lambda x: x[:-2])
# genome_seq = genome["seq"].sum()
# genes = fasta_to_dataframe(prawn_cdss_path)
# genes["name"] = genes["description"].apply(lambda x: re.findall("gene=(.*?)\]", x)[0])
# genes["protein"] = genes["description"].apply(lambda x: re.findall("protein=(.*?)\]", x)).apply(lambda x: x[0] if x else "")
# genes["location"] = genes["description"].apply(lambda x: re.findall("location=(.*?)\]", x)[0])
# genes["exons"] = genes["location"].apply(lambda x: re.findall("\d+", x)).apply(lambda x: [int(e) for e in x])
# genes["start"] = genes["exons"].apply(min)
# genes["end"] = genes["exons"].apply(max)
# genes["len"] = genes["seq"].str.len()
# genes = genes.groupby("name").apply(lambda x: x.nlargest(1, "len")).reset_index(drop=True)
# genes["chromosome_genbank_id"] = genes["id"].apply(lambda x: re.findall("N(?:C|T|W)_\d+\.\d+", x)).apply(lambda x: x[0] if x else "")
# def get_full_seq(gene):
#     chromosome = genome[genome["id"] == gene["chromosome_genbank_id"]]
#     if chromosome.empty:
#         return ""
#     sequence = chromosome["seq"].values[0]
#     return sequence[gene["start"]:gene["end"]]
# genes["full_seq"] = genes.apply(get_full_seq, axis=1)
# orf_seq = genes["seq"].sum()
# full_gene_seq = genes["full_seq"].sum()
# genome_seq = genome["seq"].sum()
def calc_competition_feats(sites,genome,genes):
    for ref_seq_label, ref_seq in {"orf": orf_seq, "full_gene": full_gene_seq, "genome": genome_seq}.items():
        for suffix_len in [5, 10, 15, 20]:
            print(f"Processing {ref_seq_label} {suffix_len}")
            sequences_with_ngg_pam = sites["sequence"].str[-suffix_len:] + "NGG"
            pool = mp.Pool(10)
            match_counts = pool.map(functools.partial(get_site_count, search_string=ref_seq), sequences_with_ngg_pam.values)
            pool.terminate()
            sites[[f"sites_count_ngg_pam_{ref_seq_label}_-{suffix_len}_d{d}" for d in range(max_mismaches)]] = pd.DataFrame(match_counts)
    
    for ref_seq_label, ref_seq in {"orf": orf_seq, "full_gene": full_gene_seq, "genome": genome_seq}.items():
        for suffix_len in [5, 10, 15, 20]:
            print(f"Processing {ref_seq_label} {suffix_len}")
            sequences = sites["sequence"].str[-suffix_len:]
            pool = mp.Pool(10)
            match_counts = pool.map(functools.partial(get_site_count, search_string=ref_seq), sequences.values)
            pool.terminate()
            sites[[f"sites_count_{ref_seq_label}_-{suffix_len}_d{d}" for d in range(max_mismaches)]] = pd.DataFrame(match_counts)
    sites[[f"sites_count_ngg_pam_genome_-20_d{d}" for d in range(max_mismaches)]].iloc[0]
    return sites