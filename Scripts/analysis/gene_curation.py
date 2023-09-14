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
import sys
import numpy as np
import pandas as pd
import pickle


# %%
import matplotlib.pyplot as plt
from plotnine import *
import seaborn as sns; sns.set_theme()
from matplotlib.colors import LogNorm

# %%
# import functions
if os.getcwd().split('/')[-1] == 'vale':
    sys.path.append('Scripts/function/') # path to use when run the script from snakemake pipeline
else: 
    sys.path.append('../function/') # path to use when run the script locally from jupyterlab

from prediction_function import *
from postprocessing_function import *
from analysis_function import *

# %%
# read in snakemake object
file = open("/s/project/vale/driver_prediction_published_202309/processed_data/snakemake/postprocessing.p",'rb')
snakemake = pickle.load(file)

# %%
experiment_path = snakemake.input.experimentDesign
features_partial_path = snakemake.input.features_partial
features_full_path = snakemake.input.features_full
res_pre_dir = os.path.dirname(snakemake.input.result)
res_post_dir = os.path.dirname(snakemake.output.result_post)

project_dir = snakemake.params.projectPath
intogen_dir = snakemake.params.intogenDir
outrider_dir = snakemake.params.outriderDir
absplice_path = snakemake.params.absplice
        
cgc_cancer_gene_path = snakemake.params.cgc_cancer_gene_processed
intogen_cancer_gene_path = snakemake.params.intogen_cancer_gene

# %%
cgc_cancer_gene = pd.read_csv(cgc_cancer_gene_path, sep='\t')
cgc_leukemia_gene = cgc_cancer_gene[["L" in x for x in cgc_cancer_gene['TissueType']]]
cgc_AML_gene = cgc_cancer_gene[["AML" in str(x) for x in cgc_cancer_gene['TumourTypes(Somatic)']]]
# cgc_cancer_gene.head()

intogen_cancer_gene = pd.read_csv(intogen_cancer_gene_path, sep='\t').rename(columns={'SYMBOL':'GeneSymbol'})
intogen_AML_gene = intogen_cancer_gene[intogen_cancer_gene.COHORT == "TCGA_WXS_LAML"]

intogen_summary = intogen_cancer_gene[['GeneSymbol', 'CANCER_TYPE', 'ROLE']].drop_duplicates()
intogen_summary = pd.merge(
    intogen_summary.groupby('GeneSymbol')['CANCER_TYPE'].apply(','.join).reset_index(),
    intogen_summary[['GeneSymbol', 'ROLE']].drop_duplicates(), 
    on='GeneSymbol')
# intogen_cancer_gene.head()

# %% [markdown]
# # define result to curate

# %%
label_gene_list = 'CGC_leukemia_gene'
sample_group = 'leukemia_14group'
model_method = 'rf'
intogen_input_feature = 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv'
outlier_input_feature = 'or,ac,absplice,fr'
coess_input_feature = ''

# %%
# read in exp
experiment_design = pd.read_csv(experiment_path, sep='\t')
experiment_design = experiment_design.fillna('')
experiment_design = add_result_paths(experiment_design, project_dir, intogen_dir)

# %%
experiment_design_group = experiment_design[
    (np.isin(experiment_design.label_gene_list, label_gene_list)) *
    (np.isin(experiment_design.model_method, model_method)) *
    (np.isin(experiment_design.intogen_input_feature, intogen_input_feature)) *
    (np.isin(experiment_design.outlier_input_feature, outlier_input_feature)) *
    (np.isin(experiment_design.coess_input_feature, coess_input_feature))
]

experiment_design_group.shape

# %%
# read in result_post
for experiment_design_group_index, experiment_design_group_row in experiment_design_group.iterrows():
    # read in res
    res_post_sample_group = pd.read_csv(experiment_design_group_row['res_post_path'], sep='\t')
    res_post_sample_group = res_post_sample_group.assign(sample_group = experiment_design_group_row['sample_group'])
    res_post_sample_group = res_post_sample_group[['GeneID', 'GeneSymbol', 'sample_group', 'Label', 'Prediction_post', 'Rank_post']]

    if 'res_post_group' in locals():
        res_post_group = pd.concat([res_post_group, res_post_sample_group])
    else:
        res_post_group = res_post_sample_group
        
# res_post_group

# %%
# read in vep
for experiment_design_group_index, experiment_design_group_row in experiment_design_group.iterrows():
    # read in res
    if os.path.exists(experiment_design_group_row['vep_path']):
        vep_sample_group = pd.read_csv(experiment_design_group_row['vep_path'], sep='\t')
        vep_sample_group[['col1', 'sampleID', 'ref', 'alt', 'coordinate']] = vep_sample_group['#Uploaded_variation'].str.split('__', expand=True)
        vep_sample_group[['chromosome', 'coordinate_range']] = vep_sample_group['Location'].str.split(':', 2, expand=True)
        vep_sample_group['coordinate']=pd.to_numeric(vep_sample_group.coordinate)
        vep_sample_group = vep_sample_group.assign(sample_group = experiment_design_group_row['sample_group'])
        vep_sample_group = vep_sample_group[['Gene', 'SYMBOL', 'chromosome', 'coordinate', 'sampleID', 'sample_group', 'IMPACT', 'Consequence', '#Uploaded_variation', 'Location']]

        if 'vep_full' in locals():
            vep_full = pd.concat([vep_full, vep_sample_group])
        else:
            vep_full = vep_sample_group
            
# vep_full

# %%
# read in dndscv
for experiment_design_group_index, experiment_design_group_row in experiment_design_group.iterrows():
    # read in res
    if os.path.exists(experiment_design_group_row['dndscv_path']):
        dndscv_sample_group = pd.read_csv(experiment_design_group_row['dndscv_path'], sep='\t')
        dndscv_sample_group = dndscv_sample_group.assign(sample_group = experiment_design_group_row['sample_group'])
        dndscv_sample_group['ROLE'] = dndscv_sample_group.apply(lambda row: set_role(row), axis=1)
        dndscv_sample_group = dndscv_sample_group[['gene_name', 'sample_group', 'ROLE']]

        if 'dndscv_full' in locals():
            dndscv_full = pd.concat([dndscv_full, dndscv_sample_group])
        else:
            dndscv_full = dndscv_sample_group
        
# dndscv_full

# %%
experiment_design_mode = experiment_design[
    (np.isin(experiment_design.label_gene_list, ['CGC_leukemia_OCG', 'CGC_leukemia_TSG'])) *
    (np.isin(experiment_design.model_method, model_method)) *
    (np.isin(experiment_design.intogen_input_feature, intogen_input_feature)) *
    (np.isin(experiment_design.outlier_input_feature, outlier_input_feature)) *
    (np.isin(experiment_design.coess_input_feature, coess_input_feature)) 
]

experiment_design_mode.shape

# %%
# read in result_mode
for experiment_design_group_index, experiment_design_group_row in experiment_design_mode.iterrows():
    # read in res
    res_post_mode_sample_group = pd.read_csv(experiment_design_group_row['res_post_path'], sep='\t')
    res_post_mode_sample_group = res_post_mode_sample_group.assign(label_gene_list = experiment_design_group_row['label_gene_list'])
    res_post_mode_sample_group = res_post_mode_sample_group.assign(sample_group = experiment_design_group_row['sample_group'])
    res_post_mode_sample_group = res_post_mode_sample_group[['GeneID', 'GeneSymbol', 'label_gene_list', 'sample_group', 'Rank_post']]

    if 'res_post_mode' in locals():
        res_post_mode = pd.concat([res_post_mode, res_post_mode_sample_group])
    else:
        res_post_mode = res_post_mode_sample_group
        
res_post_mode

# %%
experiment_design_heatmap = experiment_design[
    (np.isin(experiment_design.label_gene_list, label_gene_list)) *
    (np.isin(experiment_design.model_method, model_method)) *
    (np.isin(experiment_design.intogen_input_feature, [intogen_input_feature, ''])) *
    (np.isin(experiment_design.outlier_input_feature, [outlier_input_feature, ''])) *
    (np.isin(experiment_design.coess_input_feature, [coess_input_feature, '']))
]

experiment_design_heatmap.shape

# %%
# read in result_post for heatmap
for experiment_design_sub_index, experiment_design_sub_row in experiment_design_heatmap.iterrows():
    # read in res
    res_post_sample_group = pd.read_csv(experiment_design_sub_row['res_post_path'], sep='\t')
    res_post_sample_group = res_post_sample_group.assign(Prediction_post_scale = res_post_sample_group.Prediction_post/res_post_sample_group.Prediction_post.max())
    res_post_sample_group = res_post_sample_group.assign(sample_group = experiment_design_sub_row['sample_group'])
    res_post_sample_group = res_post_sample_group.assign(intogen_input_feature = experiment_design_sub_row['intogen_input_feature'])
    res_post_sample_group = res_post_sample_group.assign(outlier_input_feature = experiment_design_sub_row['outlier_input_feature'])
    res_post_sample_group = res_post_sample_group.assign(coess_input_feature = experiment_design_sub_row['coess_input_feature'])
    res_post_sample_group = res_post_sample_group[[
        'GeneID', 'GeneSymbol', 'Prediction_post', 'Prediction_post_scale',  'Rank_post',
        'sample_group', 'intogen_input_feature', 'outlier_input_feature', 'coess_input_feature'
    ]]
    
    if 'res_post_heatmap' in locals():
        res_post_heatmap = pd.concat([res_post_heatmap, res_post_sample_group])
    else:
        res_post_heatmap = res_post_sample_group
        
res_post_heatmap

# %%
experiment_design_spider = experiment_design[
    (np.isin(experiment_design.label_gene_list, label_gene_list)) *
    (np.isin(experiment_design.model_method, model_method)) *
    (np.isin(experiment_design.sample_group, 'leukemia_14group')) *
    (np.isin(experiment_design.intogen_input_feature, ['fml', 'clustl', 'hotmaps', 'smregions', 'cbase', 'dndscv', 'mutpanning',''])) *
    (np.isin(experiment_design.outlier_input_feature, ['or', 'ac', 'absplice', 'fr',''])) *
    (np.isin(experiment_design.coess_input_feature, ['coess_cluster','']))
]

experiment_design_spider.shape

# %%
# read in result_post for spider
for experiment_design_sub_index, experiment_design_sub_row in experiment_design_spider.iterrows():
    # read in res
    res_post_sample_group = pd.read_csv(experiment_design_sub_row['res_post_path'], sep='\t')
    res_post_sample_group = res_post_sample_group.assign(Prediction_post_scale = res_post_sample_group.Prediction_post/res_post_sample_group.Prediction_post.max())
    res_post_sample_group = res_post_sample_group.assign(sample_group = experiment_design_sub_row['sample_group'])
    res_post_sample_group = res_post_sample_group.assign(intogen_input_feature = experiment_design_sub_row['intogen_input_feature'])
    res_post_sample_group = res_post_sample_group.assign(outlier_input_feature = experiment_design_sub_row['outlier_input_feature'])
    res_post_sample_group = res_post_sample_group.assign(coess_input_feature = experiment_design_sub_row['coess_input_feature'])
    res_post_sample_group = res_post_sample_group[[
        'GeneID', 'GeneSymbol', 'Prediction_post', 'Prediction_post_scale',
        'sample_group', 'intogen_input_feature', 'outlier_input_feature', 'coess_input_feature'
    ]]    
    
    if 'res_post_spider' in locals():
        res_post_spider = pd.concat([res_post_spider, res_post_sample_group])
    else:
        res_post_spider = res_post_sample_group
        
# res_post_spider

# %%
# read in features full
features_full = pd.read_csv(features_full_path, sep='\t')

# %% [markdown]
# ### sort by affected cell type                
#     'AML', 'MDS', 'CML', 'CMML', 'Myo_group', 'Mast_group', # myeloid
#     'ALL_group', 'CLL', 'HZL_group', 'Lym_group', 'T_NHL', 'MACS_group', # lymphoid
#     'NK' # undefined

# %% [markdown]
# ### sort by tissue
#     'AML', 'MDS', 'CML', 'CMML', 'Myo_group', 'Mast_group', 'ALL_group', 'CLL', 'HZL_group', 'Lym_group', 'T_NHL', 'NK' # peripheral blood or bone marrow
#     'MACS_group' # MACS

# %% [markdown]
# # gene curation

# %% [markdown]
# ## define gene

# %%
gene_to_curate = "AFF3"

# %% [markdown]
# ## raw feature

# %%
features_curate = features_full[features_full.geneSymbol == gene_to_curate].transpose().reset_index()
features_curate.columns = ['feature', 'value']
features_curate

# %%
feature_to_check_regex = 'or-leu|ac-leu'

features_curate.loc[features_curate.feature.str.contains(feature_to_check_regex, regex=True)]

# %% [markdown]
# ## mode of action

# %%
# from dndscv - variant based
dndscv_full[dndscv_full.gene_name == gene_to_curate]

# %%
# from prediction rank post
res_post_mode_sub = res_post_mode[res_post_mode.GeneSymbol == gene_to_curate]

res_post_mode_sub = res_post_mode_sub.pivot_table(index=['GeneID', 'GeneSymbol', 'sample_group'], columns='label_gene_list',values='Rank_post').reset_index()
res_post_mode_sub.columns.name = None
res_post_mode_sub['Role'] = res_post_mode_sub.apply(lambda row: estimate_role(row), axis=1)

res_post_mode_sub

# %% [markdown]
# ## prediction post visualization

# %%
# get heatmap result
res_post_heatmap_sub = res_post_heatmap.loc[res_post_heatmap.GeneSymbol == gene_to_curate]

res_post_heatmap_sub.replace({
      #'fml,clustl,hotmaps,smregions,cbase,dndscv,mutpanning': 'genomic', 
      'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv': 'genomic',
     'or,ac,absplice,fr': 'outlier', 
     'coess_cluster': 'coess', 
     '': np.nan}, inplace = True)
res_post_heatmap_sub = res_post_heatmap_sub.assign(input_feature = res_post_heatmap_sub[['intogen_input_feature', 'outlier_input_feature', 'coess_input_feature']].apply(
    lambda x: '+'.join(x.dropna()), axis=1))

res_post_heatmap_sub

# %% [markdown]
# ### heatmap - raw prediction post

# %%
# get heatmap df
heatmap_pred_post_df = res_post_heatmap_sub.pivot_table(index=['sample_group'], columns='input_feature', values='Prediction_post')
heatmap_pred_post_df.index.name = None
heatmap_pred_post_df.columns.name = None

heatmap_pred_post_df.index = pd.CategoricalIndex(
    heatmap_pred_post_df.index, 
    categories=['leukemia_14group', 
                'AML', 'BCP_ALL', 'CLL', 'CML', 'CMML', 'HZL_group', 'Lym_group',
                'MACS_group', 'MDS', 'MDS_MPN_group', 'MPN', 'Mast_group', 'NK',
                'TL_group'
               ], ordered=True)
heatmap_pred_post_df.sort_index(level=0, inplace=True)

heatmap_pred_post_df = heatmap_pred_post_df.transpose()

heatmap_pred_post_df.index = pd.CategoricalIndex(
    heatmap_pred_post_df.index, 
    categories=['genomic', 'outlier', 'coess', 'genomic+outlier', 'genomic+coess', 'outlier+coess', 'genomic+outlier+coess'], ordered=True)

heatmap_pred_post_df.sort_index(level=0, inplace=True)

heatmap_pred_post_df

# %%
fig, ax = plt.subplots(figsize=(7,4.9))
sns.set(font_scale=1.4)
heatmap_pred_post = sns.heatmap(heatmap_pred_post_df, cmap='Reds', ax=ax, vmin=0, vmax=1)
plt.title('Raw prediction post ' + gene_to_curate, fontsize=20)
plt.show()

# %%
fig, ax = plt.subplots(figsize=(7,4.9))  
heatmap_pred_post = sns.heatmap(heatmap_pred_post_df, norm=LogNorm(), cmap='Reds', ax=ax, vmin=0, vmax=1)
plt.title('Raw prediction post ' + gene_to_curate, fontsize=20)
plt.show()

# %%
# heatmap_pred_post.get_figure().savefig('heatmap_pred_post.png', dpi=400)

# %% [markdown]
# ### heatmap - scaled prediction post

# %%
# get heatmap df
heatmap_pred_post_scaled_df = res_post_heatmap_sub.pivot_table(index=['sample_group'], columns='input_feature', values='Prediction_post_scale')
heatmap_pred_post_scaled_df.index.name = None
heatmap_pred_post_scaled_df.columns.name = None

heatmap_pred_post_scaled_df.index = pd.CategoricalIndex(
    heatmap_pred_post_scaled_df.index, 
    categories=['leukemia_14group', 
                'AML', 'BCP_ALL', 'CLL', 'CML', 'CMML', 'HZL_group', 'Lym_group',
                'MACS_group', 'MDS', 'MDS_MPN_group', 'MPN', 'Mast_group', 'NK',
                'TL_group'
               ], ordered=True)
heatmap_pred_post_scaled_df.sort_index(level=0, inplace=True)

heatmap_pred_post_scaled_df = heatmap_pred_post_scaled_df.transpose()

heatmap_pred_post_scaled_df.index = pd.CategoricalIndex(
    heatmap_pred_post_scaled_df.index, 
    categories=['genomic', 'outlier', 'coess', 'genomic+outlier', 'genomic+coess', 'outlier+coess', 'genomic+outlier+coess'], ordered=True)

heatmap_pred_post_scaled_df.sort_index(level=0, inplace=True)

# heatmap_pred_post_scaled_df

# %%
fig, ax = plt.subplots(figsize=(7,4.9))
heatmap_pred_post_scaled = sns.heatmap(heatmap_pred_post_scaled_df, cmap='Reds', ax=ax, vmin=0, vmax=1)
plt.title('Scaled prediction post ' + gene_to_curate, fontsize=20)
plt.show()

# %%
fig, ax = plt.subplots(figsize=(7,4.9))  
heatmap_pred_post_scaled = sns.heatmap(heatmap_pred_post_scaled_df, norm=LogNorm(), cmap='Reds', ax=ax, vmin=0, vmax=1)
plt.title('Scaled prediction post ' + gene_to_curate, fontsize=20)
plt.show()

# %%
# heatmap_pred_post_scaled.get_figure().savefig('heatmap_pred_post_scaled.png', dpi=400)

# %% [markdown]
# ### heatmap - rank post

# %%
# get heatmap df with ranks
heatmap_pred_post_rank = res_post_heatmap_sub.pivot_table(index=['sample_group'], columns='input_feature', values='Rank_post')
heatmap_pred_post_rank.index.name = None
heatmap_pred_post_rank.columns.name = None

heatmap_pred_post_rank.index = pd.CategoricalIndex(
    heatmap_pred_post_rank.index, 
    categories=['leukemia_14group', 
                'AML', 'BCP_ALL', 'CLL', 'CML', 'CMML', 'HZL_group', 'Lym_group',
                'MACS_group', 'MDS', 'MDS_MPN_group', 'MPN', 'Mast_group', 'NK',
                'TL_group'
               ], ordered=True)
heatmap_pred_post_rank.sort_index(level=0, inplace=True)

heatmap_pred_post_rank = heatmap_pred_post_rank.transpose()

heatmap_pred_post_rank.index = pd.CategoricalIndex(
    heatmap_pred_post_rank.index, 
    categories=['genomic', 'outlier', 'coess', 'genomic+outlier', 'genomic+coess', 'outlier+coess', 'genomic+outlier+coess'], ordered=True)

heatmap_pred_post_rank.sort_index(level=0, inplace=True)

# heatmap_pred_post_rank

# %%
fig, ax = plt.subplots(figsize=(7,4.9))  
heatmap_pred_post_rank = sns.heatmap(heatmap_pred_post_rank, cmap='Reds_r', ax=ax, annot=True, fmt='g', annot_kws={"fontsize":6})
plt.title('Rank post ' + gene_to_curate, fontsize=20)
plt.show()

# %%
# heatmap_pred_post_rank.get_figure().savefig('heatmap_pred_post_rank.png', dpi=400, bbox_inches='tight')

# %% [markdown]
# ### spider - raw prediction post

# %%
# get spider result
res_post_spider_sub = res_post_spider.loc[res_post_spider.GeneSymbol == gene_to_curate]
res_post_spider_sub.replace({'': np.nan}, inplace = True)

res_post_spider_sub = res_post_spider_sub.assign(input_feature = res_post_spider_sub[['intogen_input_feature', 'outlier_input_feature', 'coess_input_feature']].apply(
    lambda x: '+'.join(x.dropna()), axis=1))

# %%
spider_features = res_post_spider_sub.input_feature.values
spider_features = [*spider_features, spider_features[0]]

spider_pred_post = res_post_spider_sub.Prediction_post.values
spider_pred_post = [*spider_pred_post, spider_pred_post[0]]

spider_pred_post_scale = res_post_spider_sub.Prediction_post_scale.values
spider_pred_post_scale = [*spider_pred_post_scale, spider_pred_post_scale[0]]

label_loc = np.linspace(start=0, stop=2 * np.pi, num=len(spider_pred_post)) # label locations in radians

# %%
p_spider_pred_post = plt.figure(figsize=(5, 5))
plt.subplot(polar=True)
plt.ylim([0.0, 1.07])
plt.plot(label_loc, spider_pred_post)
plt.title('Single input - pred post: '+gene_to_curate, size=16)
lines, labels = plt.thetagrids(np.degrees(label_loc), labels=spider_features) # places names on label locations
plt.xticks(fontsize=14)
plt.yticks(fontsize=11)
plt.show()

# %%
# p_spider_pred_post.savefig('spider_pred_post.png', dpi=400)

# %% [markdown]
# ### spider - scaled prediction post

# %%
p_spider_pred_post_scale = plt.figure(figsize=(5, 5))
plt.subplot(polar=True)
plt.ylim([0.0, 1.07])
plt.plot(label_loc, spider_pred_post_scale)
plt.title('Single input - pred post scaled: '+gene_to_curate, size=16)
lines, labels = plt.thetagrids(np.degrees(label_loc), labels=spider_features) # places names on label locations
plt.xticks(fontsize=14)
plt.yticks(fontsize=11)
plt.show()

# %%
# p_spider_pred_post_scale.savefig('spider_pred_post_scale.png', dpi=400)

# %% [markdown]
# ## variant visualization

# %%
vep_sub = vep_full.loc[vep_full.SYMBOL == gene_to_curate]
vep_sub = vep_sub[vep_sub.SYMBOL == gene_to_curate][['sampleID', 'Gene', 'SYMBOL', 'chromosome', 'coordinate', 'Consequence', 'IMPACT', 'sample_group']]
vep_sub.IMPACT = pd.Categorical(vep_sub.IMPACT, categories=['HIGH', 'MODERATE', 'LOW'], ordered=True)
vep_sub.sample_group = pd.Categorical(
    vep_sub.sample_group, 
    categories=['leukemia_14group', 
                'AML', 'BCP_ALL', 'CLL', 'CML', 'CMML', 'HZL_group', 'Lym_group',
                'MACS_group', 'MDS', 'MDS_MPN_group', 'MPN', 'Mast_group', 'NK',
                'TL_group'
               ], ordered=True)
vep_sub

# %%
# to avoid error regarding numpy64
np.float = float  

# %%
p_vep_full = (ggplot(vep_sub, aes(x='coordinate', fill='Consequence'))
      + geom_histogram(bins = 100)
      + facet_wrap('IMPACT', ncol=1)
      
#       + scale_y_discrete(breaks = np.array([0, 5, 10, 15]))
                         
      + ggtitle(gene_to_curate)
      + xlab("Genomic position")
      + ylab("Count of variants")
      
      + theme_bw()
      # + theme(axis_title_x=element_blank(),
      + theme(axis_text_x=element_blank())
      # + theme(axis_text_y=element_text(size=12))
      + theme(legend_text = element_text(size=11))
     )

p_vep_full

# %%
p_vep_cohort = (ggplot(vep_sub, aes(x='coordinate', fill='Consequence'))
      + geom_histogram(bins = 100)
      + facet_grid('sample_group ~ IMPACT')
      
#       + scale_y_discrete(breaks = np.array([0, 5, 10, 15]))
                         
      + ggtitle(gene_to_curate)
      + xlab("Genomic position")
      + ylab("Count of variants")
      
      + theme_bw()
      + theme(axis_title_x=element_blank(),
              axis_text_x=element_blank())
     #  + theme(axis_text_y=element_text(size=12))
     )

p_vep_cohort

# %%
p_vep_cohort.save(filename = 'variant_distribution.png', height=15, width=6, units = 'in', dpi=300)

# %%
p_vep_aml = (ggplot(vep_sub.loc[vep_sub.sample_group=='AML'], aes(x='coordinate', fill='Consequence'))
      + geom_histogram(bins = 100)
      + facet_wrap('IMPACT', ncol=1)
      
#       + scale_y_discrete(breaks = np.array([0, 5, 10, 15]))
                         
      + ggtitle(gene_to_curate)
      + xlab("Genomic position")
      + ylab("Count of variants")
      
      + theme_bw()
      + theme(axis_title_x=element_blank(),
              axis_text_x=element_blank())
     #  + theme(axis_text_y=element_text(size=12))
     )

p_vep_aml

# %%
# p_vep_aml.save(filename = 'variant_distribution.png', height=3, width=4, units = 'in', dpi=300)

# %% [markdown]
# ## coess module

# %%
features_coess = get_coess_features(features_full, coess_input_features='coess_cluster')
features_coess = pd.concat([features_full[['gene_id', 'geneSymbol']], features_coess], axis=1)
features_coess = add_cancer_driver_info_feat(features_coess, cgc_cancer_gene, cgc_leukemia_gene, cgc_AML_gene, intogen_cancer_gene)
features_coess['isLabel'] = np.isin(features_coess['gene_id'], res_post_group.loc[res_post_group.Label, 'GeneID'])

# %%
coess_curate = features_curate.loc[features_curate.value!=0]
coess_curate = coess_curate.loc[coess_curate.apply(lambda x: 'coess_cluster' in x['feature'], axis=1)]
coess_curate = calculate_coess_curate(coess_curate, features_coess)
coess_curate['feature'] = coess_curate['feature'].str.replace("coess_cluster-", "")

coess_curate = coess_curate.sort_values(by=['cgc_ratio'], ascending=False)

coess_curate

# %%
coess_curate_ratio_vis = (ggplot(coess_curate, aes(x='ocg_ratio', y='tsg_ratio', size='size'))
      + geom_abline(color='grey')
      + geom_point(alpha = 0.25)
                          
      + coord_fixed(ratio=1.0)
      + xlim(0,coess_curate[['ocg_ratio', 'tsg_ratio']].max().max())
      + ylim(0,coess_curate[['ocg_ratio', 'tsg_ratio']].max().max())
                          
      + ggtitle('Co-Ess modules containing gene '+ gene_to_curate)
      + labs(size="module size")
      )

coess_curate_ratio_vis

# %%
# coess_curate_ratio_vis.save(filename = 'coess_curate_ratio_vis.png', height=3, width=3, units = 'in', dpi=300)

# %%
cluster_to_curate = 'coess_cluster-'+'d0.9-1236'

cluster_curate = features_coess.loc[features_coess[cluster_to_curate] == 1][['geneSymbol',cluster_to_curate, 'isCGC', 'isIntOGen', 'isLabel']]
cluster_curate = pd.merge(cluster_curate, cgc_cancer_gene, how='left', left_on='geneSymbol', right_on='GeneSymbol')
cluster_curate = pd.merge(cluster_curate, intogen_summary, how='left', left_on='geneSymbol', right_on='GeneSymbol')

cluster_curate[[
    'geneSymbol', 'isLabel',
    'isCGC', 'TumourTypes(Somatic)','TumourTypes(Germline)', 'CancerSyndrome', 
    'TissueType', 'RoleinCancer', 'MutationTypes', 
    'isIntOGen', 'CANCER_TYPE', 'ROLE'
]].fillna('')

# %% [markdown]
# # sample curation

# %% [markdown]
# ## define sample

# %%
sample_id_to_curate = 'MLL_15868'

# %% [markdown]
# ## varaint curation

# %%
vep_full

# %%
# var of the sample
vep_sample = vep_full[vep_full.sampleID == sample_id_to_curate]
vep_sample = vep_sample.rename(columns={'Gene':'GeneID', 'SYMBOL':'GeneSymbol'})
vep_sample.IMPACT = pd.Categorical(vep_sample.IMPACT, categories=['HIGH', 'MODERATE', 'LOW'], ordered=True)
vep_sample = vep_sample.sort_values(by=['IMPACT'])
vep_sample = add_cancer_driver_info(vep_sample, cgc_cancer_gene, cgc_leukemia_gene, cgc_AML_gene, intogen_cancer_gene)
vep_sample = vep_sample.replace(False, '')
vep_sample[0:10]

# %%
vep_sample.Consequence.value_counts()

# %%
vep_sample.loc[vep_sample.isCGC==True]

# %%
# var of the sample and the gene
vep_sample[vep_sample.GeneSymbol == gene_to_curate]
