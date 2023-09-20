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


# %%
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/extract_intogen_score.p"
pickle.dump(snakemake, open(snakemake_path, "wb"))

# %%
# file = open("/s/project/vale/driver_prediction_202303_2/processed_data/snakemake/extract_intogen_score.p",'rb')
# snakemake = pickle.load(file)

# %%
input_dir = snakemake.input.inputDir
projectPath = snakemake.params.projectPath
output_dir = projectPath + "/processed_data/intogen_feature/score"

# %%
# create output dir if not exist
if not os.path.isdir(output_dir): 
    os.mkdir(projectPath + "/processed_data/intogen_feature")
    os.mkdir(projectPath + "/processed_data/intogen_feature/score")

# %% [markdown]
# # hotmaps -log(q_value)

# %%
score_path = output_dir + "/hotmaps.tsv"

# %%
hotmaps_input_paths = glob.glob( input_dir + '/steps/hotmaps/*.out.gz')
# hotmaps_input_paths

# %%
hotmaps_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group'])

for i in hotmaps_input_paths:
#     print(i)
    hotmaps_sub_df = pd.read_csv(i, sep = '\t')
    hotmaps_sub_df = hotmaps_sub_df.loc[hotmaps_sub_df['q-value'] < 1]
 

    if hotmaps_sub_df.shape[0] != 0: 
        hotmaps_sub_df['score'] = - np.log(hotmaps_sub_df['q-value'])
        intogen_group_name = re.search('.*/MLL_WGS_MLL_(.+).out.gz', i).group(1)
        print(intogen_group_name)
        hotmaps_sub_df['group'] = intogen_group_name
        hotmaps_sub_df = hotmaps_sub_df.rename(columns = {'GENE': 'gene_name'})
        hotmaps_sub_df = hotmaps_sub_df[['gene_name', 'score', 'group']]
        
        hotmaps_score_df = pd.concat([hotmaps_score_df, hotmaps_sub_df]) 

# %%
hotmaps_score_df

# %%
hotmaps_score_df.to_csv(score_path, sep = '\t')

# %% [markdown]
# # clustl - Score sum

# %%
score_path = output_dir + "/clustl.tsv"

# %%
clustl_input_paths = glob.glob( input_dir + '/steps/oncodriveclustl/*/clusters_results.tsv')
# clustl_input_paths

# %%
clustl_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group'])

for i in clustl_input_paths:
    # print(i)
    clustl_sub_df = pd.read_csv(i, sep = '\t')
    clustl_sub_df = clustl_sub_df.loc[clustl_sub_df['P'] < 1]
 

    if clustl_sub_df.shape[0] != 0: 
        clustl_sub_df = clustl_sub_df.rename(columns = {'SCORE': 'score'}) 
        intogen_group_name = re.search('.*/MLL_WGS_MLL_(.+)/clusters_results.tsv', i).group(1)
        print(intogen_group_name)
        clustl_sub_df['group'] = intogen_group_name
        clustl_sub_df = clustl_sub_df.rename(columns = {'SYMBOL': 'gene_name'})
        clustl_sub_df = clustl_sub_df[['gene_name', 'score', 'group']]
        clustl_sub_df = pd.pivot_table(clustl_sub_df, values='score', index=['gene_name', 'group'], aggfunc=sum).reset_index()
        
        clustl_score_df = pd.concat([clustl_score_df, clustl_sub_df])

# %%
clustl_score_df

# %%
clustl_score_df.to_csv(score_path, sep = '\t')

# %% [markdown]
# # smregions U max

# %%
score_path = output_dir + "/smregions.tsv"

# %%
smregions_input_paths = glob.glob( input_dir + '/steps/smregions/*.smregions.tsv.gz')
# smregions_input_paths

# %%
smregions_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group'])

for i in smregions_input_paths:
    # print(i)
    smregions_sub_df = pd.read_csv(i, sep = '\t')
    smregions_sub_df = smregions_sub_df.loc[smregions_sub_df['Q_VALUE'] < 1]
 

    if smregions_sub_df.shape[0] != 0: 
        smregions_sub_df = smregions_sub_df.rename(columns = {'U': 'score'})
        intogen_group_name = re.search('.*/MLL_WGS_MLL_(.+).smregions.tsv.gz', i).group(1)
        print(intogen_group_name)
        smregions_sub_df['group'] = intogen_group_name
        smregions_sub_df = smregions_sub_df.rename(columns = {'HUGO_SYMBOL': 'gene_name'})
        smregions_sub_df = smregions_sub_df[['gene_name', 'score', 'group']]
        smregions_sub_df = pd.pivot_table(smregions_sub_df, values='score', index=['gene_name', 'group'], aggfunc=max).reset_index()
        
        smregions_score_df = pd.concat([smregions_score_df, smregions_sub_df])

# %%
smregions_score_df

# %%
smregions_score_df.to_csv(score_path, sep = '\t')

# %% [markdown]
# # fml -log(q_value)

# %%
score_path = output_dir + "/fml.tsv"

# %%
fml_input_paths = glob.glob( input_dir + '/steps/oncodrivefml/out/*-oncodrivefml.tsv.gz')
fml_input_paths

# %%
fml_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group'])

for i in fml_input_paths:
    # print(i)
    fml_sub_df = pd.read_csv(i, sep = '\t')
    fml_sub_df = fml_sub_df.loc[fml_sub_df['Q_VALUE'] < 1]

    if fml_sub_df.shape[0] != 0: 
        fml_sub_df['score'] = - np.log(fml_sub_df['Q_VALUE'])
        intogen_group_name = re.search('.*/MLL_WGS_MLL_(.+)-oncodrivefml.tsv.gz', i).group(1)
        print(intogen_group_name)
        fml_sub_df['group'] = intogen_group_name
        fml_sub_df = fml_sub_df.rename(columns = {'SYMBOL': 'gene_name'})
        fml_sub_df = fml_sub_df[['gene_name', 'score', 'group']]
        
        fml_score_df = pd.concat([fml_score_df, fml_sub_df])

# %%
fml_score_df

# %%
fml_score_df.to_csv(score_path, sep = '\t')

# %% [markdown]
# # mutpanning -log(FDR) with filter .query('Significance < 0.1').query('FDR < 1')

# %%
score_path = output_dir + "/mutpanning.tsv"

# %%
mutpanning_input_paths = glob.glob( input_dir + '/steps/mutpanning/out/SignificanceFiltered/Significance*.txt')
mutpanning_input_paths

# %%
mutpanning_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group'])

for i in mutpanning_input_paths:
    # print(i)
    mutpanning_sub_df = pd.read_csv(i, sep = '\t')
    mutpanning_sub_df = mutpanning_sub_df.loc[mutpanning_sub_df['FDR'] < 1].loc[mutpanning_sub_df['Significance'] < 0.1]

    if mutpanning_sub_df.shape[0] != 0: 
        mutpanning_sub_df.loc[mutpanning_sub_df.FDR == 0, 'FDR'] = min(mutpanning_sub_df.loc[mutpanning_sub_df.FDR != 0, 'FDR']) # fill FDR 0 to the minimun FDR
        mutpanning_sub_df['score'] = - np.log(mutpanning_sub_df['FDR'])
        intogen_group_name = re.search('.SignificanceMLL_WGS_MLL_(.+).txt', i).group(1)
        print(intogen_group_name)
        mutpanning_sub_df['group'] = intogen_group_name
        mutpanning_sub_df = mutpanning_sub_df.rename(columns = {'Name': 'gene_name'})
        mutpanning_sub_df = mutpanning_sub_df[['gene_name', 'score', 'group']]
        
        mutpanning_score_df = pd.concat([mutpanning_score_df, mutpanning_sub_df])

# %%
mutpanning_score_df

# %%
mutpanning_score_df.to_csv(score_path, sep = '\t')

# %% [markdown]
# # dndscv -log(qallsubs_cv)

# %%
score_path = output_dir + "/dndscv.tsv"

# %%
dndscv_input_paths = glob.glob( input_dir + '/steps/dndscv/*.dndscv.tsv.gz')
# dndscv_input_paths

# %%
dndscv_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group'])

for i in dndscv_input_paths:
    # print(i)
    dndscv_sub_df = pd.read_csv(i, sep = '\t')

    if dndscv_sub_df.shape[0] != 0: 
        dndscv_sub_df.loc[dndscv_sub_df.qallsubs_cv == 0, 'qallsubs_cv'] = min(dndscv_sub_df.loc[dndscv_sub_df.qallsubs_cv != 0, 'qallsubs_cv']) # fill FDR 0 to the minimun FDR
        dndscv_sub_df['score'] = - np.log(dndscv_sub_df['qallsubs_cv'])
        intogen_group_name = re.search('.MLL_WGS_MLL_(.+).dndscv.tsv.gz', i).group(1)
        print(intogen_group_name)
        dndscv_sub_df['group'] = intogen_group_name
        dndscv_sub_df = dndscv_sub_df[['gene_name', 'score', 'group']]
        
        dndscv_score_df = pd.concat([dndscv_score_df, dndscv_sub_df])

# %%
dndscv_score_df

# %%
dndscv_score_df.to_csv(score_path, sep = '\t')

# %% [markdown]
# # cbase -log(q_value)

# %%
score_path = output_dir + "/cbase.tsv"

# %%
cbase_input_paths = glob.glob( input_dir + '/steps/cbase/*.cbase.tsv.gz')
# cbase_input_paths

# %%
cbase_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group'])

for i in cbase_input_paths:
    # print(i)
    cbase_sub_df = pd.read_csv(i, sep = '\t')
    cbase_sub_df = cbase_sub_df.loc[cbase_sub_df['q_pos'] < 1]
 

    if cbase_sub_df.shape[0] != 0: 
        cbase_sub_df.loc[cbase_sub_df.q_pos == 0, 'q_pos'] = min(cbase_sub_df.loc[cbase_sub_df.q_pos != 0, 'q_pos']) # fill FDR 0 to the minimun FDR
        cbase_sub_df['score'] = - np.log(cbase_sub_df['q_pos'])
        intogen_group_name = re.search('.*/MLL_WGS_MLL_(.+).cbase.tsv.gz', i).group(1)
        print(intogen_group_name)
        cbase_sub_df['group'] = intogen_group_name
        cbase_sub_df = cbase_sub_df.rename(columns = {'gene': 'gene_name'})
        cbase_sub_df = cbase_sub_df[['gene_name', 'score', 'group']]
        
        cbase_score_df = pd.concat([cbase_score_df, cbase_sub_df])

# %%
cbase_score_df

# %%
cbase_score_df.to_csv(score_path, sep = '\t')
