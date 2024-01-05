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
# import modules
import os
import re
import glob
import pickle
import pandas as pd
import numpy as np


# %%
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/intogen_tab.p"
pickle.dump(snakemake, open(snakemake_path, "wb"))

# %%
# file = open("/s/project/vale/driver_prediction_202304/processed_data/snakemake/intogen_tab.p",'rb')
# snakemake = pickle.load(file)

# %%
input_dir = snakemake.input.inputDir
out_path = snakemake.output.intogen_tab

# %%
# subset single group
manuscript_wording = pd.read_csv(snakemake.params.manuscriptWording, sep='\t')
manuscript_wording.columns = [i.replace(" ", "_") for i in manuscript_wording.columns]

manuscript_wording['Study_group_filename'] = manuscript_wording['Study_group']
manuscript_wording['Study_group_filename'] = [i.replace(" ", "_") for i in manuscript_wording.Study_group_filename]
manuscript_wording['Study_group_filename'] = [i.replace("-", "_") for i in manuscript_wording.Study_group_filename]
manuscript_wording['Study_group_filename'] = [i.replace("/", "_") for i in manuscript_wording.Study_group_filename]
 
manuscript_wording['Study_group_during_analysis_upper'] = manuscript_wording['Study_group_during_analysis'].str.upper()
manuscript_wording

# %% [markdown]
# # hotmaps

# %%
hotmaps_input_paths = glob.glob( input_dir + '/steps/hotmaps/*.out.gz')
hotmaps_input_paths

# %%
hotmaps_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group', 'metric'])

for i in hotmaps_input_paths:
#     print(i)
    hotmaps_sub_df = pd.read_csv(i, sep = '\t') 

    if hotmaps_sub_df.shape[0] != 0: 
        hotmaps_sub_df['score'] = hotmaps_sub_df['q-value']
        intogen_group_name = re.search('.*/MLL_WGS_MLL_(.+).out.gz', i).group(1)
        print(intogen_group_name)
        hotmaps_sub_df['group'] = intogen_group_name
        hotmaps_sub_df['metric'] = 'q_value'
        hotmaps_sub_df = hotmaps_sub_df.rename(columns = {'GENE': 'gene_name'})
        hotmaps_sub_df = hotmaps_sub_df[['gene_name', 'score', 'group', 'metric']]
        
        hotmaps_score_df = pd.concat([hotmaps_score_df, hotmaps_sub_df]) 

# %%
hotmaps_score_df

# %% [markdown]
# # clustl

# %%
clustl_input_paths = glob.glob( input_dir + '/steps/oncodriveclustl/*/clusters_results.tsv')
clustl_input_paths

# %%
clustl_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group', 'metric'])

for i in clustl_input_paths:
    # print(i)
    clustl_sub_df = pd.read_csv(i, sep = '\t')

    if clustl_sub_df.shape[0] != 0: 
        clustl_sub_df = clustl_sub_df.rename(columns = {'SCORE': 'score'}) 
        intogen_group_name = re.search('.*/MLL_WGS_MLL_(.+)/clusters_results.tsv', i).group(1)
        print(intogen_group_name)
        clustl_sub_df['group'] = intogen_group_name
        clustl_sub_df['metric'] = 'SCORE'
        clustl_sub_df = clustl_sub_df.rename(columns = {'SYMBOL': 'gene_name'})
        clustl_sub_df = clustl_sub_df[['gene_name', 'score', 'group', 'metric']]
        clustl_sub_df = pd.pivot_table(clustl_sub_df, values='score', index=['gene_name', 'group', 'metric'], aggfunc=sum).reset_index()
        
        clustl_score_df = pd.concat([clustl_score_df, clustl_sub_df])

# %%
clustl_score_df

# %% [markdown]
# # smregions

# %%
smregions_input_paths = glob.glob( input_dir + '/steps/smregions/*.smregions.tsv.gz')
smregions_input_paths

# %%
smregions_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group', 'metric'])

for i in smregions_input_paths:
    # print(i)
    smregions_sub_df = pd.read_csv(i, sep = '\t')

    if smregions_sub_df.shape[0] != 0: 
        smregions_sub_df = smregions_sub_df.rename(columns = {'U': 'score'})
        intogen_group_name = re.search('.*/MLL_WGS_MLL_(.+).smregions.tsv.gz', i).group(1)
        print(intogen_group_name)
        smregions_sub_df['group'] = intogen_group_name
        smregions_sub_df['metric'] = 'U'
        smregions_sub_df = smregions_sub_df.rename(columns = {'HUGO_SYMBOL': 'gene_name'})
        smregions_sub_df = smregions_sub_df[['gene_name', 'score', 'group', 'metric']]
        smregions_sub_df = pd.pivot_table(smregions_sub_df, values='score', index=['gene_name', 'group', 'metric'], aggfunc=max).reset_index()
        
        smregions_score_df = pd.concat([smregions_score_df, smregions_sub_df])

# %%
smregions_score_df

# %% [markdown]
# # fml

# %%
fml_input_paths = glob.glob( input_dir + '/steps/oncodrivefml/out/*-oncodrivefml.tsv.gz')
fml_input_paths

# %%
fml_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group', 'metric'])

for i in fml_input_paths:
    # print(i)
    fml_sub_df = pd.read_csv(i, sep = '\t')

    if fml_sub_df.shape[0] != 0: 
        fml_sub_df['score'] = fml_sub_df['Q_VALUE']
        intogen_group_name = re.search('.*/MLL_WGS_MLL_(.+)-oncodrivefml.tsv.gz', i).group(1)
        print(intogen_group_name)
        fml_sub_df['group'] = intogen_group_name
        fml_sub_df['metric'] = 'Q_VALUE'
        fml_sub_df = fml_sub_df.rename(columns = {'SYMBOL': 'gene_name'})
        fml_sub_df = fml_sub_df[['gene_name', 'score', 'group', 'metric']]
        
        fml_score_df = pd.concat([fml_score_df, fml_sub_df])

# %%
fml_score_df

# %% [markdown]
# # mutpanning 

# %%
mutpanning_input_paths = glob.glob( input_dir + '/steps/mutpanning/out/SignificanceFiltered/Significance*.txt')
mutpanning_input_paths

# %%
mutpanning_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group', 'metric'])

for i in mutpanning_input_paths:
    # print(i)
    mutpanning_sub_df = pd.read_csv(i, sep = '\t')

    if mutpanning_sub_df.shape[0] != 0: 
        mutpanning_sub_df['score'] = mutpanning_sub_df['FDR']
        intogen_group_name = re.search('.SignificanceMLL_WGS_MLL_(.+).txt', i).group(1)
        print(intogen_group_name)
        mutpanning_sub_df['group'] = intogen_group_name
        mutpanning_sub_df['metric'] = 'FDR'
        mutpanning_sub_df = mutpanning_sub_df.rename(columns = {'Name': 'gene_name'})
        mutpanning_sub_df = mutpanning_sub_df[['gene_name', 'score', 'group', 'metric']]
        
        mutpanning_score_df = pd.concat([mutpanning_score_df, mutpanning_sub_df])

# %%
mutpanning_score_df

# %% [markdown]
# # dndscv

# %%
dndscv_input_paths = glob.glob( input_dir + '/steps/dndscv/*.dndscv.tsv.gz')
dndscv_input_paths

# %%
dndscv_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group', 'metric'])

for i in dndscv_input_paths:
    # print(i)
    dndscv_sub_df = pd.read_csv(i, sep = '\t')

    if dndscv_sub_df.shape[0] != 0: 
        dndscv_sub_df['score'] = dndscv_sub_df['qallsubs_cv']
        intogen_group_name = re.search('.MLL_WGS_MLL_(.+).dndscv.tsv.gz', i).group(1)
        print(intogen_group_name)
        dndscv_sub_df['group'] = intogen_group_name
        dndscv_sub_df['metric'] = 'qallsubs_cv'
        dndscv_sub_df = dndscv_sub_df[['gene_name', 'score', 'group', 'metric']]
        
        dndscv_score_df = pd.concat([dndscv_score_df, dndscv_sub_df])

# %%
dndscv_score_df

# %% [markdown]
# # cbase

# %%
cbase_input_paths = glob.glob( input_dir + '/steps/cbase/*.cbase.tsv.gz')
cbase_input_paths

# %%
cbase_score_df = pd.DataFrame(columns=['gene_name', 'score', 'group', 'metric'])

for i in cbase_input_paths:
    # print(i)
    cbase_sub_df = pd.read_csv(i, sep = '\t')

    if cbase_sub_df.shape[0] != 0: 
        cbase_sub_df['score'] = cbase_sub_df['q_pos']
        intogen_group_name = re.search('.*/MLL_WGS_MLL_(.+).cbase.tsv.gz', i).group(1)
        print(intogen_group_name)
        cbase_sub_df['group'] = intogen_group_name
        cbase_sub_df['metric'] = 'q_pos'
        cbase_sub_df = cbase_sub_df.rename(columns = {'gene': 'gene_name'})
        cbase_sub_df = cbase_sub_df[['gene_name', 'score', 'group', 'metric']]
        
        cbase_score_df = pd.concat([cbase_score_df, cbase_sub_df])

# %%
cbase_score_df

# %% [markdown]
# # merge

# %%
hotmaps_score_df['method']='hotmaps'
clustl_score_df['method']='clustl'
smregions_score_df['method']='smregions'
fml_score_df['method']='fml'
mutpanning_score_df['method']='mutpanning'
dndscv_score_df['method']='dndscv'
cbase_score_df['method']='cbase'

# %%
score_df = pd.concat([hotmaps_score_df,
                      clustl_score_df,
                      smregions_score_df,
                      fml_score_df,
                      mutpanning_score_df,
                      dndscv_score_df,
                      cbase_score_df])

# %%
score_df = pd.merge(score_df, manuscript_wording[['Study_group_during_analysis_upper', 'Study_group_filename']].drop_duplicates(), 
                    how = 'left', left_on = 'group', right_on = 'Study_group_during_analysis_upper')
score_df['Colname'] = score_df.apply(lambda x: '-'.join([x['Study_group_filename'], x['method'], x['metric']]), axis=1)

score_df = score_df.pivot_table(index=['gene_name'], columns=['Colname'], values = 'score', aggfunc='mean')

# %%
score_df 

# %%
score_df.to_csv(out_path)

# %%
