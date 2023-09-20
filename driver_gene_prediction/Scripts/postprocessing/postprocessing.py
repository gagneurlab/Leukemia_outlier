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
import sys
import numpy as np
import pandas as pd
import pickle


# %%
import matplotlib.pyplot as plt
from plotnine import *


# %%
# import functions
if os.getcwd().split('/')[-1] == 'vale':
    sys.path.append('Scripts/function/') # path to use when run the script from snakemake pipeline
else: 
    sys.path.append('../function/') # path to use when run the script locally from jupyterlab

from prediction_function import *
from postprocessing_function import *

# %%
# save snakemake object
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/postprocessing.p"
pickle.dump(snakemake, open(snakemake_path, "wb"))

# %%
# file = open("/s/project/vale/driver_prediction_202207/processed_data/snakemake/postprocessing.p",'rb')
# snakemake = pickle.load(file)

# %%
experiment_no = int(snakemake.params.sampID)
experiment_path = snakemake.input.experimentDesign
features_partial_path = snakemake.input.features_partial
features_full_path = snakemake.input.features_full
res_path = snakemake.input.result
res_dir = os.path.dirname(res_path)

res_post_path = snakemake.output.result_post
rp_path = snakemake.output.rank_proportion

cgc_cancer_gene_path = snakemake.params.cgc_cancer_gene_processed
intogen_cancer_gene_path = snakemake.params.intogen_cancer_gene

# %%
# experiment_no = 1931

# %% tags=[]
cgc_cancer_gene = pd.read_csv(cgc_cancer_gene_path, sep='\t')
cgc_leukemia_gene = cgc_cancer_gene[["L" in x for x in cgc_cancer_gene['TissueType']]]
cgc_AML_gene = cgc_cancer_gene[["AML" in str(x) for x in cgc_cancer_gene['TumourTypes(Somatic)']]]
# cgc_cancer_gene.head()

intogen_cancer_gene = pd.read_csv(intogen_cancer_gene_path, sep='\t').rename(columns={'SYMBOL':'GeneSymbol'})
# intogen_cancer_gene.head()

# %%
experiment_df = pd.read_csv(experiment_path, sep='\t')
experiment_df = experiment_df.fillna('')

_, sample_group, _, _ = get_param_training(experiment_no, experiment_df)
intogen_input_feature, outlier_input_feature, coess_input_feature = get_param_feature(experiment_no, experiment_df)

# %%
experiment_df.loc[experiment_df.experiment_no == experiment_no]

# %%
# extracting features 
# read in features full
features_full = pd.read_csv(features_partial_path, sep='\t')
features_post = get_features_post(features_full, sample_group, intogen_input_feature, outlier_input_feature)

# features_post

# %%
# read in res
res = pd.read_csv(res_path, sep='\t', index_col=0)
res_info = aggregate_res_by_gene(res)
res_info = add_cancer_driver_info(res_info, cgc_cancer_gene, cgc_leukemia_gene, cgc_AML_gene, intogen_cancer_gene)

# res_info

# %%
# remove coess only predictions: set prediction of pure coess gene to 0
res_info = add_prediction_rank_post(res_info, intogen_input_feature, outlier_input_feature, coess_input_feature, features_post)

# res_info

# %%
# adding prediction result of each component
res_info = add_intogen_res(res_info, intogen_input_feature, experiment_no, experiment_df, res_dir)
res_info = add_outrider_res(res_info, outlier_input_feature, experiment_no, experiment_df, res_dir)
res_info = add_coess_res(res_info, coess_input_feature, experiment_no, experiment_df, res_dir)
res_info = add_wo_intogen_res(res_info, intogen_input_feature, experiment_no, experiment_df, res_dir)
res_info = add_wo_outrider_res(res_info, outlier_input_feature, experiment_no, experiment_df, res_dir)
res_info = add_wo_coess_res(res_info, coess_input_feature, experiment_no, experiment_df, res_dir)

# res_info

# %%
# order columns
res_info = res_info[[
    'GeneID','GeneSymbol','Label',
    'Prediction','Rank', 'Prediction_post','Rank_post', 'only_coess', 
    'intogenPred', 'intogenPred_Rank', 'outlierPred', 'outlierPred_Rank', 'coessPred', 'coessPred_Rank', 
    'woIntogenPred', 'woIntogenPred_Rank', 'woOutlierPred', 'woOutlierPred_Rank', 'woCoessPred', 'woCoessPred_Rank',
    'isCGC','isLeukemia','isAML','isOCG','isTSG','isIntOGen',
    'RoleCGC','RoleIntogen'
]]

# res_info

# %%
# merge features
res_post = pd.merge(res_info, features_post, on='GeneID', how='left')

res_post

# %%
res_post.to_csv(res_post_path, sep='\t', index=False)

# %%
# calculate proportion
rp_df = calculate_rank_proportion(res_info)

rp_df

# %%
rp_df.to_csv(rp_path, sep='\t', index=False)

# %%
# res_post.to_csv('/s/project/vale/driver_prediction_202207/processed_results/post/result_2247.tsv', sep='\t')
