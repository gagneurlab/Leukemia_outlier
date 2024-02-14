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
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/figure_4_prep_rp.p"
pickle.dump(snakemake, open(snakemake_path, "wb")) 

# %%
# file = open("/s/project/vale/driver_prediction_202401/processed_data/snakemake/figure_4_prep_rp.p",'rb')
# snakemake = pickle.load(file)

# %%
# import functions
if os.getcwd().split('/')[-1] == 'driver_gene_prediction':
    sys.path.append('Scripts/function/') # path to use when run the script from snakemake pipeline
else: 
    sys.path.append('../function/') # path to use when run the script locally from jupyterlab

from analysis_function import *
from postprocessing_function import *


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

# %%
cgc_cancer_gene = pd.read_csv(snakemake.params.cgc_cancer_gene_processed, sep='\t')
cgc_leukemia_gene = cgc_cancer_gene[["L" in x for x in cgc_cancer_gene['TissueType']]]
cgc_AML_gene = cgc_cancer_gene[["AML" in str(x) for x in cgc_cancer_gene['TumourTypes(Somatic)']]]

intogen_cancer_gene = pd.read_csv(snakemake.params.intogen_cancer_gene, sep='\t').rename(columns={'SYMBOL':'GeneSymbol'})

# %% [markdown]
# # Benchmark - Precision Recall Curve/Average Precision

# %%
res_benchmark = pd.read_csv(snakemake.input.mll_benchmark, sep='\t')
res_benchmark = res_benchmark.assign(Method_Sample_group = res_benchmark.apply(lambda x:'-'.join([x['Method'], x['Sample_group']]), axis=1))
res_benchmark

# %%
res_benchmark = res_benchmark.loc[np.isin(res_benchmark.Method, 
        ['intogen', 'cbase', 'clustl', 'dndscv', 'fml', 'hotmaps', 'mutpanning','smregions'])]
res_benchmark = res_benchmark.loc[np.isin(res_benchmark.Sample_group, 
        ['AML', 'MDS', 'Lym_group'])]

res_benchmark

# %%
# calculate Rank and Proportion
rp_df = pd.DataFrame(columns=['Rank', 'Criteria', 'Proportion'])
  
for method_group in res_benchmark.Method_Sample_group.unique():
    # read in res
    res = res_benchmark[res_benchmark.Method_Sample_group == method_group]
    res = add_cancer_driver_info(res, cgc_cancer_gene, cgc_leukemia_gene, cgc_AML_gene, intogen_cancer_gene)
    res = res.sort_values(by='Prediction', ascending=False)
    res = res.assign(Rank = res.Prediction.rank(ascending=False, method='first'))
    res = res.iloc[0:200]

    # create proportion df for rank
    rp_sub_df = pd.DataFrame(columns=['isLabel', 'isCGC', 'isLeukemia', 'isOCG', 'isTSG', 'isIntOGen'], index=res.Rank.unique())

    # calculate proportion for rank
    for i in res.Rank.unique():
        rp_sub_df.loc[i, 'isLabel'] = float(res[res.Rank <= i][['Label']].sum()/i)
        rp_sub_df.loc[i, 'isCGC'] = float(res[res.Rank <= i][['isCGC']].sum()/i)
        rp_sub_df.loc[i, 'isLeukemia'] = float(res[res.Rank <= i][['isLeukemia']].sum()/i)
        rp_sub_df.loc[i, 'isOCG'] = float(res[res.Rank <= i][['isOCG']].sum()/i)
        rp_sub_df.loc[i, 'isTSG'] = float(res[res.Rank <= i][['isTSG']].sum()/i)
        rp_sub_df.loc[i, 'isIntOGen'] = float(res[res.Rank <= i][['isIntOGen']].sum()/i)

        
    # melt proportion for plotting for rank
    rp_sub_df = rp_sub_df.reset_index().rename(columns={'index':'Rank'})
    rp_sub_df['GeneID'] = res['GeneID'].tolist()
    rp_sub_df['GeneSymbol'] = res['GeneSymbol'].tolist()
    rp_sub_df = pd.melt(rp_sub_df, id_vars=['Rank', 'GeneID', 'GeneSymbol'], var_name='Criteria', value_name='Proportion')
    rp_sub_df['Criteria'] = pd.Categorical(rp_sub_df['Criteria'], 
                                           categories=['isLabel', 'isCGC', 'isLeukemia', 'isOCG', 'isTSG', 'isIntOGen'], 
                                           ordered=True)
    
    rp_sub_df['Method_Sample_group'] = method_group

    rp_df = pd.concat([rp_df, rp_sub_df])

# %%
rp_df[['Method', 'Sample_group']] = rp_df.Method_Sample_group.str.split("-", expand = True)
rp_df

# %%
rp_df.to_csv(snakemake.output.mll_benchmark_rp, sep='\t', index=False)
