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
import datetime
import numpy as np
import pandas as pd
import itertools

# %%
# this path should be identical to the file `wbuild.yaml` variable `projectPath`  
project_path =  "/s/project/vale/driver_prediction_202312" 
if not os.path.exists(project_path):
    os.mkdir(project_path)


# %%
def add_exp(experiment_df,
            label_gene_list, sample_group, model_method, 
            intogen_input_feature, outlier_input_feature, coess_input_feature, 
            random_seed, solvers, penalty, max_iter, 
            n_estimators, min_samples_split, max_depth, 
            weight_true, weight_false, tree_method):

    iterables = [ 
        label_gene_list, sample_group, 
        intogen_input_feature, outlier_input_feature, coess_input_feature,
        model_method
                ]

    for t in itertools.product(*iterables):

        experiment_df_sub = pd.DataFrame(data={

            'label_gene_list': [t[0]],
            'sample_group': [t[1]],

            'intogen_input_feature': [t[2]],
            'outlier_input_feature': [t[3]],
            'coess_input_feature': [t[4]],

            'model_method': [t[5]],
            'random_seeds': [random_seed],

            'solver': [solvers],
            'penalty': [penalty], 
            'max_iter': [max_iter],

            'n_estimators': [n_estimators],
            'min_samples_split': [min_samples_split],
            'max_depth': [max_depth],

            'weight_true': [weight_true],
            'weight_false': [weight_false],

            'tree_method': [tree_method]
        })

        experiment_df = pd.concat([experiment_df, experiment_df_sub])
        
    return(experiment_df)


# %%
# directory paths
experiment_filepath = project_path + "/experiment_design.tsv" 
experiment_dir = os.path.dirname(experiment_filepath)
experiment_filename = "experiment_df_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + ".tsv"

experiment_filepath_time = os.path.join(experiment_dir, experiment_filename)

# %%
experiment_df = pd.DataFrame(columns=[
    'label_gene_list', 'sample_group',
    'intogen_input_feature', 'outlier_input_feature', 'coess_input_feature', 
    'random_seeds', 'model_method', 
    'solver', 'penalty', 'max_iter', 
    'n_estimators', 'min_samples_split', 'max_depth',
    'weight_true', 'weight_false'
])

# %%
# default params
random_seed = 10

# logistic regression (lr)
# solvers = ['lbfgs', 'saga']
# penalty = ['none', 'l1']
solvers = 'lbfgs'
penalty = 'none'
max_iter = 100000

# random forest (rf)
n_estimators = 100
min_samples_split = 19
max_depth = 10

weight_true = 1
weight_false = 1

tree_method = 'hist'

# %%
label_gene_list = [
#     "CGC_AML_OCG",
#     "CGC_AML_TSG",
    "CGC_leukemia_OCG",
    "CGC_leukemia_TSG",
    "CGC_leukemia_gene",
    # "CGC_cancer_OCG",
    # "CGC_cancer_TSG",
    # "CGC_cancer_gene",
#     "IntOGen_cancer_OCG",
#     "IntOGen_cancer_TSG",
#     "IntOGen_cancer_gene",
    "MLL_leukemia_gene",
    "MLL_CGC_leukemia_gene"
           ]

sample_group = [
    "leukemia_14group", 
    "AML", "MDS", "Lym_group", "MPN", "CLL", 
    "MACS_group", "MDS_MPN_group", "BCP_ALL", "TL_group", "CMML", 
    "HZL_group", "Mast_group", "CML", "NK"   
]

model_method = [
    "lr", "rf", 
    # "xgb", "xgb_op"
]

# %% [markdown]
# # add intogen setup 

# %%
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
    "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv"
]
outlier_input_feature = [""]
coess_input_feature = [""]

# %%
experiment_df = add_exp(experiment_df, 
                        label_gene_list, sample_group, model_method, 
                        intogen_input_feature, outlier_input_feature, coess_input_feature, 
                        random_seed, solvers, penalty, max_iter, 
                        n_estimators, min_samples_split, max_depth, 
                        weight_true, weight_false, tree_method)

# %% [markdown]
# # add outlier setup 

# %%
intogen_input_feature = [""]
outlier_input_feature = ["or", "ac", "absplice", "fr", 
                         "or,ac", "or,absplice", "ac,absplice", "or,fr", "ac,fr", "absplice,fr",
                         "or,ac,absplice", "or,ac,fr", "ac,absplice,fr", "or,absplice,fr", 
                         "or,ac,absplice,fr"]
coess_input_feature = [""]

# %%
experiment_df = add_exp(experiment_df, 
                        label_gene_list, sample_group, model_method, 
                        intogen_input_feature, outlier_input_feature, coess_input_feature, 
                        random_seed, solvers, penalty, max_iter, 
                        n_estimators, min_samples_split, max_depth, 
                        weight_true, weight_false, tree_method)

# %% [markdown]
# # add combined setup

# %%
intogen_input_feature = ["clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv"]
outlier_input_feature = ["or", "ac", "absplice", "fr", 
                         "or,ac", "or,absplice", "ac,absplice", "or,fr", "ac,fr", "absplice,fr",
                         "or,ac,absplice", "or,ac,fr", "ac,absplice,fr", "or,absplice,fr", 
                         "or,ac,absplice,fr"]
coess_input_feature = [""]

# %%
experiment_df = add_exp(experiment_df, 
                        label_gene_list, sample_group, model_method, 
                        intogen_input_feature, outlier_input_feature, coess_input_feature, 
                        random_seed, solvers, penalty, max_iter, 
                        n_estimators, min_samples_split, max_depth, 
                        weight_true, weight_false, tree_method)

# %% [markdown] toc-hr-collapsed=true
# # add intogen leave 1 out setup 

# %%
intogen_input_feature = [
    'hotmaps,smregions,fml,cbase,mutpanning,dndscv',
    'clustl,smregions,fml,cbase,mutpanning,dndscv',
    'clustl,hotmaps,fml,cbase,mutpanning,dndscv',
    'clustl,hotmaps,smregions,cbase,mutpanning,dndscv',
    'clustl,hotmaps,smregions,fml,mutpanning,dndscv',
    'clustl,hotmaps,smregions,fml,cbase,dndscv',
    'clustl,hotmaps,smregions,fml,cbase,mutpanning'
]
outlier_input_feature = ["or,ac,absplice,fr"]
coess_input_feature = [""]

# %%
# experiment_df = add_exp(experiment_df, 
#                         label_gene_list, sample_group, model_method, 
#                         intogen_input_feature, outlier_input_feature, coess_input_feature, 
#                         random_seed, solvers, penalty, max_iter, 
#                         n_estimators, min_samples_split, max_depth, 
#                         weight_true, weight_false, tree_method)

# %% [markdown]
# # processing

# %%
# remove empty rows
empty_row = (
    (experiment_df.intogen_input_feature == "") *
    (experiment_df.outlier_input_feature == "") *
    (experiment_df.coess_input_feature == "") 
)

experiment_df = experiment_df[np.logical_not(empty_row)]

# %%
# add experiment_no
experiment_df = experiment_df.reset_index().drop('index', axis=1).reset_index().rename(columns={'index':'experiment_no'})

# %%
# check duplicates
sum(experiment_df.drop('experiment_no', axis=1).duplicated())

# %%
experiment_df.to_csv(experiment_filepath_time, sep='\t', index = False)

# %%
if os.path.exists(experiment_filepath):
    os.remove(experiment_filepath)
    
os.symlink(experiment_filepath_time, experiment_filepath)
