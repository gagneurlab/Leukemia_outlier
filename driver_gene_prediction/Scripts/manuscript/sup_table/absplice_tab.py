# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:anaconda-vale_202204]
#     language: python
#     name: conda-env-anaconda-vale_202204-py
# ---

# %%
import os
import pickle
import pandas as pd
import numpy as np
import re

from functools import reduce
from itertools import compress
from tqdm import tqdm

# %%
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/absplice_tab.p"
pickle.dump(snakemake, open(snakemake_path, "wb"))

# %%
# file = open("/s/project/vale/driver_prediction_202304/processed_data/snakemake/absplice_tab.p",'rb')
# snakemake = pickle.load(file)

# %%
gencode_path = snakemake.params.gencode
out_path = snakemake.output.absplice_tab
single_group = snakemake.params.single_group[0]
sa_path = snakemake.params.samp_anno

# %%
# subset pr coding
gencode = pd.read_csv(gencode_path)
gencode = gencode[gencode.gene_type=='protein_coding']
gencode = gencode.rename(columns={'gene_id':'gene_id_unique'})
gencode['gene_id'] = [i.split('.', 1)[0] for i in gencode.gene_id_unique]
absplice_tab = gencode[["gene_id"]]

# %%
# subset single group
samp_anno = pd.read_csv(sa_path, sep='\t')
samp_anno = samp_anno.loc[samp_anno.DROP_GROUP.str.contains(single_group)]
samp_anno['Diag'] = [i.replace(" ", "_") for i in samp_anno.Diag]
samp_anno['Diag'] = [i.replace("-", "_") for i in samp_anno.Diag]
samp_anno['Diag'] = [i.replace("/", "_") for i in samp_anno.Diag]
samp_anno.shape

# %%
absplice_full_df = pd.read_csv(snakemake.input.absplice_res, sep='\t')
absplice_full_df

# %% [markdown]
# # generate table

# %%
diag_list = samp_anno.Diag.unique()
diag_list

# %%
for i in diag_list:
    print(i)
    
    diag_id = samp_anno.loc[samp_anno.Diag==i, 'ArrayID']
    absplice_diag_df = absplice_full_df.loc[absplice_full_df.sampleID.isin(diag_id)]

    absplice_sub_df = absplice_diag_df.loc[absplice_diag_df.AbSplice_DNA >= 0.2]
    max_02_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.max, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-ab-max_0.2'})
    mean_02_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.mean, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-ab-mean_0.2'})
    freq_02_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=len, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-ab-freq_0.2'})

    absplice_sub_df = absplice_diag_df.loc[absplice_diag_df.AbSplice_DNA >= 0.05]
    max_005_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.max, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-ab-max_0.05'})
    mean_005_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.mean, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-ab-mean_0.05'})
    freq_005_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=len, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-ab-freq_0.05'})

    absplice_sub_df = absplice_diag_df
    max_001_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.max, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-ab-max_0.01'})
    mean_001_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.mean, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-ab-mean_0.01'})
    freq_001_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=len, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-ab-freq_0.01'})

    df_list = [
        max_02_df, mean_02_df, freq_02_df, 
        max_005_df, mean_005_df, freq_005_df, 
        max_001_df, mean_001_df, freq_001_df]

    merge_df = reduce(lambda  left,right: pd.merge(left,right,on=['gene_id'], how='outer'), df_list).fillna(0)
    absplice_tab = pd.merge(absplice_tab, merge_df, on='gene_id', how='left')

# %%
absplice_tab

# %%
absplice_tab.to_csv(out_path, index=False)

# %%
