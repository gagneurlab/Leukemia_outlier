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
#     display_name: Python [conda env:absplice]
#     language: python
#     name: absplice
# ---

# %% [markdown]
# # filter criteria
# - pr coding gene
# - cutoff (can be defined in code chunk 6)
# - tissue by cohort_group
# - merge by gene

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
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/read_in_absplice_res.p"
pickle.dump(snakemake, open(snakemake_path, "wb"))


# %%
# file = open("/s/project/vale/driver_prediction_published_202309/processed_data/snakemake/read_in_absplice_res.p",'rb')
# snakemake = pickle.load(file)

# %%
def get_abs_max_rows(df, groupby, max_col, dropna=True):
    return df.reset_index() \
        .sort_values(by=max_col, key=abs, ascending=False) \
        .drop_duplicates(subset=groupby) \
        .set_index(groupby)


# %% [markdown]
# # preprocess

# %%
absplice_files = snakemake.input.absplice
gencode_path = snakemake.params.gencode
out_path = snakemake.output.absplice_res
out_var_path = snakemake.output.absplice_res_var
single_group = snakemake.params.single_group[0]
sa_path = snakemake.params.samp_anno 

# %%
gencode = pd.read_csv(gencode_path)
gencode = gencode.rename(columns={'gene_id':'gene_id_unique'})
gencode['gene_id'] = [i.split('.', 1)[0] for i in gencode.gene_id_unique]
feature_absplice = gencode[["gene_id"]]

# %%
samp_anno = pd.read_csv(sa_path, sep='\t')
samp_anno_explode = samp_anno.assign(DROP_GROUP=samp_anno.DROP_GROUP.str.split(',')).explode('DROP_GROUP')

# %% [markdown]
# # read in

# %%
# subset for study group
samp_anno = samp_anno.loc[samp_anno.DROP_GROUP.str.contains(single_group)]
samp_anno.shape

# %%
# subset for existing 
samp_anno = samp_anno[np.logical_not(samp_anno.DNA_ABSPLICE_FILE.isna())]
samp_anno.shape

# %%
cutoff = 0.01

# %%
absplice_var_df = pd.DataFrame(columns=[
    'gene_id', 'gene_name', 'sampleID', 'variant', 'tissue', 'delta_logit_psi',
    'delta_psi', 'delta_score', 'splice_site_is_expressed', 'AbSplice_DNA',
    'gene_type'])

for absplice_path in tqdm(samp_anno.DNA_ABSPLICE_FILE): 
    # read in

    df = pd.read_parquet(absplice_path)
    df = df.reset_index()

    # subset for protein coding
    df = pd.merge(df, gencode[['gene_id', 'gene_type', 'gene_name']], on='gene_id')
    df = df.loc[df.gene_type=='protein_coding']

    # get sampleID
    sampID = samp_anno.loc[samp_anno.DNA_ABSPLICE_FILE==absplice_path].ArrayID.values[0]
    cohort_group = samp_anno.loc[samp_anno.DNA_ABSPLICE_FILE==absplice_path].Cohort_group.values[0]

    # subset for tissue
    df = df.loc[df.tissue == cohort_group]
    df = df.assign(sampleID = sampID)



    # subset by cutoff
    df_sub = df[
        (df['AbSplice_DNA'] >= cutoff)
    ].sort_values(by='AbSplice_DNA', ascending=False)
    df_sub = df_sub.reset_index()
    
    absplice_var_df = pd.concat([absplice_var_df, df_sub])

# %%
absplice_var_df

# %%
absplice_var_df.shape

# %%
absplice_var_df.to_csv(out_var_path, sep='\t', index=False)

# %%
# aggregate by gene
absplice_gene_df = get_abs_max_rows(
    df = absplice_var_df.set_index(['gene_id', 'sampleID']),
    groupby = ['gene_id', 'sampleID'],
    max_col = 'AbSplice_DNA',
#     dropna=False
)
absplice_gene_df = absplice_gene_df.reset_index()

# %%
absplice_gene_df

# %%
absplice_gene_df.shape

# %%
absplice_gene_df.to_csv(out_path, sep='\t', index=False)

# %%
