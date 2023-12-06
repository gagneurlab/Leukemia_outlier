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
import time
import pickle
import numpy as np
import pandas as pd


# %%
# import statistic modules
from sklearn.model_selection import StratifiedKFold


# %%
# import functions
if os.getcwd().split('/')[-1] == 'driver_gene_prediction':
    sys.path.append('Scripts/function/') # path to use when run the script from snakemake pipeline
else: 
    sys.path.append('../function/') # path to use when run the script locally from jupyterlab
    
from prediction_function import *

# %%
# save snakemake object
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/prediction.p"
pickle.dump(snakemake, open(snakemake_path, "wb"))

# %%
# file = open("/s/project/vale/driver_prediction_202303_2/processed_data/snakemake/prediction.p",'rb')
# snakemake = pickle.load(file)

# %%
experiment_no = int(snakemake.params.sampID)
label_dir_path = os.path.dirname(snakemake.input.label[1])
experiment_path = snakemake.input.experimentDesign
features_partial_path = snakemake.input.features_partial
features_full_path = snakemake.input.features_full
res_full_path = snakemake.output.result
coeff_full_path = snakemake.output.coeff

# %%
# experiment_no = 621

# %%
experiment_df = pd.read_csv(experiment_path, sep='\t')
experiment_df = experiment_df.fillna('')

label_gene_list, sample_group, random_seeds, model_method = get_param_training(experiment_no, experiment_df)
intogen_input_feature, outlier_input_feature, coess_input_feature = get_param_feature(experiment_no, experiment_df)
solver, penalty, max_iter = get_param_lr(experiment_no, experiment_df)
n_estimators, min_samples_split, max_depth, weight_true, weight_false = get_param_rf(experiment_no, experiment_df)
tree_method = get_param_xgb(experiment_no, experiment_df)

# %%
experiment_df.loc[experiment_df.experiment_no == experiment_no]

# %%
# read in features full
if coess_input_feature == '':
    features_full = pd.read_csv(features_partial_path, sep='\t')
else:
    features_full = pd.read_csv(features_full_path, sep='\t')

# %%
# generate features_full and features_model
features_model = get_feature_model(features_full, sample_group, intogen_input_feature, outlier_input_feature, coess_input_feature)

features_model

# %% tags=[]
# intialize result table
print_res_path(res_full_path)
res_df = initialize_res_df()

# intialize coefficiant table
print_coeff_path(coeff_full_path)
coeff_df = intitialize_coeff_df(model_method)

# read in label gene list
driver_genes = read_in_label_gene(label_dir_path, label_gene_list)

# create labels
labels = create_label(features_full, driver_genes)

# %%
# start training
start_time = start_training()

for random_round in range(random_seeds):

    start_round(random_round)
    
    # initialize prediction array
    predictions = intialize_prediction_array(features_full)

    # Using Skicit-learn for cross validation
    skf = StratifiedKFold(n_splits = 5, shuffle = True, random_state = random_round)
    
    for train_index_num, test_index_num in skf.split(features_full, labels):

        if model_method == "lr":
            predictions_test, coeff_df_sub = run_logistic_regression(features_model, labels,
                                                                     train_index_num, test_index_num,
                                                                     solver, penalty, max_iter, random_round)
            
        if model_method == "rf":
            predictions_test, coeff_df_sub = run_random_forest(features_model, labels,
                                                               train_index_num, test_index_num,
                                                               n_estimators, min_samples_split, max_depth,
                                                               weight_true, weight_false, random_round)
            
        if model_method == "xgb":
            predictions_test, coeff_df_sub = run_xgboost(features_model, labels,
                                                         train_index_num, test_index_num,
                                                         tree_method, random_round)
            
        if model_method == "xgb_op":
            predictions_test, coeff_df_sub = run_xgboost_optimized(features_model, labels,
                                                                   train_index_num, test_index_num,
                                                                   tree_method, random_round)
            
        # concat predictioins and coeff_df
        predictions[test_index_num] = predictions_test
        coeff_df_sub['random_repeat'] = random_round
        coeff_df = pd.concat([coeff_df, coeff_df_sub]) 

    # create res_df
    res_df_sub = generate_result_table(features_full, labels, predictions, random_round)
    res_df = pd.concat([res_df, res_df_sub])

# %%
# save result
coeff_df.to_csv(coeff_full_path, sep='\t', index=False)    
res_df.to_csv(res_full_path, sep='\t', index=False)

# %%
end_training(start_time)

# %%
