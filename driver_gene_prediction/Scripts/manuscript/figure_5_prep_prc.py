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
import pickle
import itertools
import numpy as np
import pandas as pd
from functools import reduce


# %%
from plotnine import *
import matplotlib.pyplot as plt


# %%
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve

# %%
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/figure_5_prep_prc.p"
pickle.dump(snakemake, open(snakemake_path, "wb")) 

# %%
# file = open("/s/project/vale/driver_prediction_202304/processed_data/snakemake/figure_5_prep_prc.p",'rb')
# snakemake = pickle.load(file)

# %%
# import functions
if os.getcwd().split('/')[-1] == 'driver_gene_prediction':
    sys.path.append('Scripts/function/') # path to use when run the script from snakemake pipeline
else: 
    sys.path.append('../function/') # path to use when run the script locally from jupyterlab

from analysis_function import *


# %%
def add_input_feature(exp_viz):

    exp_viz_sub = exp_viz[['intogen_input_feature', 'outlier_input_feature', 'coess_input_feature']]
    exp_viz = exp_viz.assign(input_feature = exp_viz_sub.apply(lambda x: "-".join(x.replace("", np.nan).dropna()), axis=1))
    
    return exp_viz


# %% [markdown]
# # Define plots

# %%
experiment_path = snakemake.params.experimentDesign
project_dir = snakemake.params.projectPath
intogen_dir = snakemake.params.intogenDir
single_group = snakemake.params.single_group 

# %%
label_gene_list = 'CGC_leukemia_gene'
sample_group = single_group
model_method = 'rf'
intogen_input_feature = [
    "clustl",
    "hotmaps",
    "smregions",
    "fml",
    "cbase",
    "mutpanning",
    "dndscv",
    "clustl,hotmaps",
    "clustl,hotmaps,smregions",
    "clustl,hotmaps,smregions,fml",
    "clustl,hotmaps,smregions,fml,cbase",
    "clustl,hotmaps,smregions,fml,cbase,mutpanning",
    "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv", 
    ""]
outlier_input_feature = ["or", "ac", "absplice", "fr", 
                         "or,ac", "or,absplice", "ac,absplice", "or,fr", "ac,fr", "absplice,fr",
                         "or,ac,absplice", "or,ac,fr", "ac,absplice,fr", "or,absplice,fr", 
                         "or,ac,absplice,fr", 
                         ""]
coess_input_feature = ['']

# %% [markdown]
# # Preparation

# %%
experiment_design = pd.read_csv(experiment_path, sep='\t')
experiment_design = experiment_design.fillna("")

# %%
exp_viz = experiment_design[
    (np.isin(experiment_design.label_gene_list, label_gene_list)) *
    (np.isin(experiment_design.sample_group, sample_group)) *
    (np.isin(experiment_design.model_method, model_method)) *
    (np.isin(experiment_design.intogen_input_feature, intogen_input_feature)) *
    (np.isin(experiment_design.outlier_input_feature, outlier_input_feature)) *
    (np.isin(experiment_design.coess_input_feature, coess_input_feature))
]
exp_viz = add_input_feature(exp_viz)
exp_viz = add_result_paths(exp_viz, project_dir, intogen_dir)

exp_viz

# %% [markdown]
# # MLL - Precision Recall Curve

# %%
# calculate Recall and Precision
performance_df = pd.DataFrame(columns=['Training_setup', 'random_repeat', 'AP']) 
prc_df = pd.DataFrame(columns=['Training_setup', 'random_repeat', 'precision','recall'])
prc_ci_df = pd.DataFrame(columns=['Training_setup', 'recall', 'precision_mean', 'precision_median', 'precision_std'])
    
for exp_viz_index, exp_viz_row in exp_viz.iterrows():
    # read in res
    res = pd.read_csv(exp_viz_row['res_pre_path'], sep='\t')
    input_feature = exp_viz_row['input_feature']
    
    performance_df_ts = pd.DataFrame(columns=['Training_setup', 'random_repeat', 'AP']) 

    for rep in res.random_repeat.unique():
        
        # extract res of certain random_repeat
        res_sub = res[res['random_repeat'] == rep]
        labels = res_sub['Label']
        predictions = res_sub['Prediction']

        # calculate AP
        average_precision = average_precision_score(labels, predictions)
        precision, recall, _ = precision_recall_curve(labels, predictions)

        performance_df_rep = pd.DataFrame({'Training_setup': input_feature, 
                                           'random_repeat': rep,
                                           'AP': average_precision}, index=[0])
        performance_df_ts = pd.concat([performance_df_ts, performance_df_rep])
        performance_df = pd.concat([performance_df, performance_df_rep])
        
        prc_df_rep = pd.DataFrame({'Training_setup': input_feature, 
                                   'random_repeat': rep,
                                   'precision': precision,
                                   'recall': recall}) # unique for each random repeat
        prc_df = pd.concat([prc_df, prc_df_rep]) # unique for each training setup

    prc_df_ts = prc_df[prc_df.Training_setup == input_feature]
    for rc in prc_df_ts.recall.unique():
        
        # using all the precision per recall
        pr_per_rc = prc_df_ts[prc_df_ts.recall == rc].precision
        prc_ci_df_rc = pd.DataFrame({'Training_setup': input_feature, 
                                     'recall': rc, 
                                     'precision_mean': pr_per_rc.mean(),
                                     'precision_median': pr_per_rc.median(),
                                     'precision_std': pr_per_rc.quantile(),
                                     'precision_up': pr_per_rc.quantile(0.9),
                                     'precision_dn': pr_per_rc.quantile(0.1)}, index=[0])
        prc_ci_df = pd.concat([prc_ci_df, prc_ci_df_rc]) 

# %%
prc_ci_df

# %%
prc_ci_df.Training_setup.value_counts()

# %%
prc_ci_df.to_csv(snakemake.output.mll_prc, sep='\t', index=False)

# %%
performance_df.to_csv(snakemake.output.mll_ap_full, sep='\t', index=False)

# %% [markdown]
# # MLL - Average Precision plot

# %%
# calculate median, 90% q, 10% q performance
ap_dfs =  [pd.pivot_table(performance_df, values='AP', index = ['Training_setup'], aggfunc=np.median).reset_index().rename(columns={'AP':'AP_median'}),
           pd.pivot_table(performance_df, values='AP', index = ['Training_setup'], aggfunc=lambda y: np.percentile(y, 90)).reset_index().rename(columns={'AP':'AP_9q'}),
           pd.pivot_table(performance_df, values='AP', index = ['Training_setup'], aggfunc=lambda y: np.percentile(y, 10)).reset_index().rename(columns={'AP':'AP_1q'})]
ap_df = reduce(lambda  left,right: pd.merge(left,right,on=['Training_setup']), ap_dfs)

# %%
ap_df

# %%
ap_df.to_csv(snakemake.output.mll_ap, sep='\t', index=False)

# %% [markdown]
# # Benchmark - Precision Recall Curve/Average Precision

# %%
res_benchmark = pd.read_csv(snakemake.input.mll_benchmark, sep='\t')
res_benchmark = res_benchmark.assign(Method_Sample_group = res_benchmark.apply(lambda x:'-'.join([x['Method'], x['Sample_group']]), axis=1))
res_benchmark

# %%
# calculate Recall and Precision
performance_df = pd.DataFrame(columns=['Method_Sample_group', 'AP']) 
prc_df = pd.DataFrame(columns=['Method_Sample_group', 'precision','recall'])
  
for method_group in res_benchmark.Method_Sample_group.unique():
    # read in res
    res = res_benchmark[res_benchmark.Method_Sample_group == method_group]
    labels = res['Label']
    predictions = res['Prediction']

    # calculate AP
    average_precision = average_precision_score(labels, predictions)
    precision, recall, _ = precision_recall_curve(labels, predictions)

    performance_df_ts = pd.DataFrame({'Method_Sample_group': method_group,
                                      'AP': average_precision}, index=[0])
    

    prc_df_rep = pd.DataFrame({'Method_Sample_group': method_group, 
                               'precision': precision,
                               'recall': recall}) # unique for each random repeat
    
    prc_df = pd.concat([prc_df, prc_df_rep]) # unique for each training setup
    performance_df = pd.concat([performance_df, performance_df_ts])

# %%
prc_df[['Method', 'Sample_group']] = prc_df.Method_Sample_group.str.split("-", expand = True)
prc_df

# %%
prc_df.to_csv(snakemake.output.mll_benchmark_prc, sep='\t', index=False)

# %% [markdown]
# # Benchmark - Average Precision plot

# %%
performance_df[['Method', 'Sample_group']] = performance_df.Method_Sample_group.str.split("-", expand = True)
performance_df

# %%
performance_df.to_csv(snakemake.output.mll_benchmark_ap, sep='\t', index=False)

# %%

# %%
