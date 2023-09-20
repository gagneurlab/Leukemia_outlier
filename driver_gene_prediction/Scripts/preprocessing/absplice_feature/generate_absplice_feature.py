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
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/generate_absplice_feature.p"
pickle.dump(snakemake, open(snakemake_path, "wb"))

# %%
# file = open("/s/project/vale/driver_prediction_202304/processed_data/snakemake/generate_absplice_feature.p",'rb')
# snakemake = pickle.load(file)

# %% [markdown]
# # preprocess

# %%
gencode_path = snakemake.params.gencode
out_path = snakemake.output.absplice
samp_groups = snakemake.params.samp_groups
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
absplice_full_df = pd.read_csv(snakemake.input.absplice_res, sep='\t')
absplice_full_df

# %% [markdown]
# # generate feature

# %%
# for groups
for i in samp_groups:
    print(i)
    
    sample_group_id = samp_anno_explode.loc[samp_anno_explode.DROP_GROUP==i, 'ArrayID'].unique()
    absplice_group_df = absplice_full_df.loc[absplice_full_df.sampleID.isin(sample_group_id)]

    absplice_sub_df = absplice_group_df.loc[absplice_group_df.AbSplice_DNA >= 0.2]
    max_02_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.max, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-max_02'})
    mean_02_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.mean, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-mean_02'})
    freq_02_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=len, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-freq_02'})

    absplice_sub_df = absplice_group_df.loc[absplice_group_df.AbSplice_DNA >= 0.05]
    max_005_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.max, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-max_005'})
    mean_005_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.mean, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-mean_005'})
    freq_005_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=len, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-freq_005'})

    absplice_sub_df = absplice_group_df
    max_001_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.max, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-max_001'})
    mean_001_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.mean, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-mean_001'})
    freq_001_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=len, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': f'{i}-freq_001'})

    df_list = [
        max_02_df, mean_02_df, freq_02_df, 
        max_005_df, mean_005_df, freq_005_df, 
        max_001_df, mean_001_df, freq_001_df]

    merge_df = reduce(lambda  left,right: pd.merge(left,right,on=['gene_id'], how='outer'), df_list).fillna(0)
    feature_absplice = pd.merge(feature_absplice, merge_df, on='gene_id', how='left')

# %%
# for the single group
single_group

absplice_sub_df = absplice_full_df.loc[absplice_full_df.AbSplice_DNA >= 0.2]
max_02_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.max, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': single_group + '-max_02'})
mean_02_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.mean, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': single_group + '-mean_02'})
freq_02_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=len, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': single_group + '-freq_02'})

absplice_sub_df = absplice_full_df.loc[absplice_full_df.AbSplice_DNA >= 0.05]
max_005_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.max, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': single_group + '-max_005'})
mean_005_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.mean, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': single_group + '-mean_005'})
freq_005_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=len, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': single_group + '-freq_005'})

absplice_sub_df = absplice_full_df
max_001_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.max, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': single_group + '-max_001'})
mean_001_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=np.mean, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': single_group + '-mean_001'})
freq_001_df = pd.pivot_table(absplice_sub_df, values='AbSplice_DNA', index=['gene_id'], aggfunc=len, fill_value=0).reset_index().rename(columns={'AbSplice_DNA': single_group + '-freq_001'})

df_list = [
    max_02_df, mean_02_df, freq_02_df, 
    max_005_df, mean_005_df, freq_005_df, 
    max_001_df, mean_001_df, freq_001_df]

merge_df = reduce(lambda  left,right: pd.merge(left,right,on=['gene_id'], how='outer'), df_list).fillna(0)
feature_absplice = pd.merge(feature_absplice, merge_df, on='gene_id', how='left')

# %%
feature_absplice.columns = feature_absplice.columns[:1].tolist() + [f'absplice-{col}' for col in feature_absplice.columns[1:]]
feature_absplice = feature_absplice.fillna(0)

feature_absplice

# %%
feature_absplice.to_csv(out_path, sep='\t', index=False)
