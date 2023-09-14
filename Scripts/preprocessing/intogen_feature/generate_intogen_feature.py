# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python [conda env:anaconda-vale_202204]
#     language: python
#     name: conda-env-anaconda-vale_202204-py
# ---

# %%
# import modules
import os
import re
import glob
import pickle
import pandas as pd
import numpy as np
import math


# %%
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/generate_intogen_feature.p"
pickle.dump(snakemake, open(snakemake_path, "wb"))

# %%
# file = open("/s/project/vale/driver_prediction_202303_2/processed_data/snakemake/generate_intogen_feature.p",'rb')
# snakemake = pickle.load(file)

# %%
intogen_dir = snakemake.params.intogenDir
gencode_path = snakemake.params.gencode
single_group = snakemake.params.single_group[0]
input_dir = os.path.dirname(snakemake.input.score)
sa_path = snakemake.params.samp_anno

# %%
gencode = pd.read_csv(gencode_path).rename(columns = {'gene_id':'gene_id_unique'})
gencode['gene_id'] = gencode['gene_id_unique'].str.split(pat=".", expand=True)[0]

# %% tags=[]
samp_anno = pd.read_csv(sa_path, sep='\t')
samp_anno_explode = samp_anno.assign(DROP_GROUP=samp_anno.DROP_GROUP.str.split(',')).explode('DROP_GROUP')
samp_anno_explode['DROP_GROUP_CAP'] = samp_anno_explode['DROP_GROUP'].apply(str.upper)

# %%
output_dir = input_dir.replace('score', 'feature')

# create output dir if not exist
if not os.path.isdir(output_dir): 
    os.mkdir(output_dir)


# %% [markdown]
# # preparation

# %%
def set_role(data, distance_threshold=0.1):
    """Set the role according to the DNDS output"""
    if data['wmis_cv'] < 1 and data['wnon_cv'] < 1:  # threshold
        return "ambiguous"
    
    # Check wmis
    wmis = data['wmis_cv']
    if wmis >= 1 and data["n_mis"] == 0:
        wmis = 1

    # Check wnon
    wnon = data['wnon_cv']
    if wnon >= 1 and data["n_non"] == 0:
        wnon = 1
        
    # Those cases with w_non and w_mis <=1 are not informative
    if wnon <= 1 and wmis <= 1:
        return "ambiguous"

    distance = (wmis - wnon) / math.sqrt(2)
    if distance_threshold is not None and abs(distance) < distance_threshold:
        return "ambiguous"
    else:
        if distance > 0:
            return 'Act'
        elif distance < 0:
            return 'LoF'
        else:
            return "ambiguous"


# %%
def feature_combine_by_cohort_intogen(feature_df, single_group_name):
    feature_output_df = feature_df[['gene_id']]
    feature_output_df = feature_output_df.assign(score_up = feature_df.filter(like='_up', axis=1).sum(axis=1))
    feature_output_df = feature_output_df.assign(score_dn = feature_df.filter(like='_dn', axis=1).sum(axis=1))
    feature_output_df = feature_output_df.assign(score_am = feature_df.filter(like='_am', axis=1).sum(axis=1))

    
    feature_output_df.columns = feature_output_df.columns[:1].tolist() + [single_group_name + '-' + str(col) for col in feature_output_df.columns[1:]]
    return feature_output_df

# %% [markdown]
# # read in data

# %%
score_paths = glob.glob( input_dir + '/*.tsv')
score_paths

# %%
for score_path in score_paths:
    
    output_path = score_path.replace('score', 'feature')
    
    score_input_df = pd.read_csv(score_path, sep = '\t', index_col=0)
    score_input_df = pd.merge(score_input_df, 
                              samp_anno_explode[['DROP_GROUP_CAP', 'DROP_GROUP']].drop_duplicates(), 
                              left_on='group', right_on='DROP_GROUP_CAP', how='left')
    score_input_df = score_input_df.drop(['group', 'DROP_GROUP_CAP'], axis=1).rename(columns={'DROP_GROUP':'group'})
    sample_groups = score_input_df.group.unique()

    for sample_group in sample_groups:
        print(sample_group)

        score_input_sub_df = score_input_df.loc[score_input_df.group == sample_group]

        # get role df
        dndscv_path = intogen_dir + '/steps/dndscv/MLL_WGS_MLL_' + sample_group.upper() + '.dndscv.tsv.gz'
        dndscv_sub_df = pd.read_csv(dndscv_path, sep = '\t')
        dndscv_sub_df['ROLE'] = dndscv_sub_df.apply(lambda row: set_role(row), axis=1)

        # add role
        score_input_sub_df = score_input_sub_df.merge(dndscv_sub_df[["gene_name", "ROLE"]], on='gene_name', how='left')
        
        # add gene_id
        score_input_sub_df = score_input_sub_df.merge(gencode[["gene_name", "gene_id"]], on='gene_name', how='left')
        
        # remove na
        if sum(score_input_sub_df.gene_id.isna()) != 0:  
            print('gene_name not found in gencode:')
            print(score_input_sub_df.gene_name[score_input_sub_df.gene_id.isna()])
            score_input_sub_df = score_input_sub_df.loc[np.logical_not(score_input_sub_df.gene_id.isna())]
        
        # pivot table
        score_input_sub_df = pd.pivot_table(score_input_sub_df, values = ['score'], index=['gene_id'], columns=['ROLE'], aggfunc=sum).reset_index().fillna(0)
#         print(score_input_sub_df.head())
        score_input_sub_df.columns = ['_'.join(col) for col in score_input_sub_df.columns.values]
        score_input_sub_df.columns = score_input_sub_df.columns[:1].tolist() + [sample_group + '-' + str(col) for col in score_input_sub_df.columns[1:]]
        score_input_sub_df = score_input_sub_df.rename(columns = {'gene_id_':'gene_id'})
#         print(score_input_sub_df.head()) 
            
        # merge    
        if not ('score_output_df' in locals() or 'score_output_df' in globals()):
            score_output_df = score_input_sub_df
        else:
            score_output_df = score_output_df.merge(score_input_sub_df, on='gene_id', how='outer')

    # fix na and column name        
    score_output_df = score_output_df.fillna(0)
    score_output_df.columns = score_output_df.columns.str.replace(
        'score_Act', 'score_up', regex=False).str.replace(
        'score_LoF', 'score_dn', regex=False).str.replace(
        'score_ambiguous', 'score_am', regex=False)
    
    # create combined column
    score_combind_df = feature_combine_by_cohort_intogen(score_output_df, single_group) # line to combine features of different cohort
    score_output_df = pd.merge(score_output_df, score_combind_df, on = 'gene_id')
    
    # rename column
    method = re.search('./score/(.+).tsv', score_path).group(1)
    score_output_df.columns = score_output_df.columns[:1].tolist() + [method + '-' + str(col) for col in score_output_df.columns[1:]]

    score_output_df.to_csv(output_path, sep = '\t', index=False)
    print(output_path)
    print(score_output_df.head()) 
    
    del score_output_df

# %%

# %%
