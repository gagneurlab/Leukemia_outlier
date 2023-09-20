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
import glob
import pickle
import pandas as pd
import numpy as np


# %%
from functools import reduce

# %%
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/merge_feature.p"
pickle.dump(snakemake, open(snakemake_path, "wb"))

# %%
file = open("/s/project/vale/driver_prediction_published_202309/processed_data/snakemake/merge_feature.p",'rb')
snakemake = pickle.load(file)

# %%
gencode_path = snakemake.params.gencode
intogen_dir = os.path.dirname(snakemake.input.intogen)
outlier_dir = snakemake.params.outlierDir
absplice_dir = os.path.dirname(snakemake.input.absplice)
coess_modules_dir = os.path.dirname(snakemake.input.coess_modules)
joint_embedding_1_path = snakemake.input.joint_embedding_1
features_partial_path = snakemake.output.features_partial
features_full_path = snakemake.output.features_full

# %%
# gencode_path = "/s/project/mll/drop_2020apr/processed_data/aberrant_expression/v33b/gene_name_mapping_v33b.tsv"
# intogen_dir = '/s/project/vale/driver_prediction_202207/processed_data/intogen_feature/feature'
# outlier_dir = '/s/project/vale/driver_prediction_202207/processed_data/outlier_feature'
# coess_modules_dir = '/s/project/vale/driver_prediction_202207/processed_data/coess_feature'
# output_path = '/s/project/vale/driver_prediction_202207/processed_data/features_full.tsv'

# %%
# create empty features 
gencode = pd.read_csv(gencode_path)
gencode = gencode[gencode.gene_type == 'protein_coding'] # subset to protein coding
gencode = gencode.rename(columns={'gene_id':'gene_id_unique'})
gencode['gene_id'] = [i.split('.', 1)[0] for i in gencode.gene_id_unique]
features = gencode[['gene_id', 'gene_name_orig']].drop_duplicates().rename(columns={'gene_name_orig':'geneSymbol'})
# features

# %%
# intOGen
print("")
print("intOGen feature dimensions:")

features_intogen_paths = glob.glob(intogen_dir + '/*.tsv')
    
for features_intogen_path in features_intogen_paths:
    features_intogen_sub = pd.read_csv(features_intogen_path, sep="\t")
    print(features_intogen_sub.shape)
    
    if not ('features_intogen' in locals() or 'features_intogen' in globals()):
        features_intogen = features_intogen_sub
    else:
        features_intogen = features_intogen.merge(features_intogen_sub, on='gene_id', how='outer')
        
# features_intogen

# %%
(features_intogen.fillna(0) != 0).any(axis=0).value_counts()-1

# %%
# abSplice
print("")
print("abSplice feature dimensions:")

features_absplice_paths = glob.glob(absplice_dir + '/*.tsv')
    
for features_absplice_path in features_absplice_paths:
    features_absplice_sub = pd.read_csv(features_absplice_path, sep="\t")
    print(features_absplice_sub.shape)
    
    if not ('features_absplice' in locals() or 'features_absplice' in globals()):
        features_absplice = features_absplice_sub
    else:
        features_absplice = features_absplice.merge(features_absplice_sub, on='gene_id', how='outer')
        
# features_absplice

# %%
(features_absplice.fillna(0) != 0).any(axis=0).value_counts()-1

# %%
# outlier
print("")
print("outlier feature dimensions:")

features_outlier_paths = glob.glob(outlier_dir + '/*.tsv')
    
for features_outlier_path in features_outlier_paths:
    features_outlier_sub = pd.read_csv(features_outlier_path, sep="\t")
    print(features_outlier_sub.shape)
    
    if not ('features_outlier' in locals() or 'features_outlier' in globals()):
        features_outlier = features_outlier_sub
    else:
        features_outlier = features_outlier.merge(features_outlier_sub, on='gene_id', how='outer')
        
# features_outlier

# %%
outlier_dir

# %%
features_or = features_outlier.loc[:, features_outlier.columns.str.contains('or-', regex=True)]
(features_or.fillna(0) != 0).any(axis=0).value_counts()

# %%
features_ac = features_outlier.loc[:, features_outlier.columns.str.contains('ac-', regex=True)]
(features_ac.fillna(0) != 0).any(axis=0).value_counts()

# %%
features_fr = features_outlier.loc[:, features_outlier.columns.str.contains('fr-', regex=True)]
(features_fr.fillna(0) != 0).any(axis=0).value_counts()

# %%
# Co-essentiality modules
print("")
print("Co-essentiality module feature dimensions:")

features_coess_modules_paths = glob.glob(coess_modules_dir + '/coess_cluster*.tsv')
    
for features_coess_modules_path in features_coess_modules_paths:
    features_coess_modules_sub = pd.read_csv(features_coess_modules_path, sep="\t")
    print(features_coess_modules_sub.shape)
    
    if not ('features_coess_modules' in locals() or 'features_coess_modules' in globals()):
        features_coess_modules = features_coess_modules_sub
    else:
        features_coess_modules = features_coess_modules.merge(features_coess_modules_sub, on='gene_id', how='outer')
        
# features_coess_modules

# %%
(features_coess_modules.fillna(0) != 0).any(axis=0).value_counts()-1

# %%
274+126+308+161+308+5224

# %%
# Co-essentiality embedding 1
print("")
print("Co-essentiality embedding_1 feature dimensions:")
    
features_joint_embedding_1 = pd.read_csv(joint_embedding_1_path, sep="\t")
print(features_joint_embedding_1.shape)

# features_joint_embedding_1

# %%
# rename columns
features_partial_dic = {
    "feature": features, 
    "intogen": features_intogen,
    "absplice": features_absplice,
    "outlier": features_outlier
    }        

# %%
features_partial_list = list(features_partial_dic.values())
features_partial = reduce(lambda  left,right: pd.merge(left,right,on='gene_id', how='left'), features_partial_list).fillna(0)
features_partial = features_partial.loc[:, (features_partial != 0).any(axis=0)]
features_partial = features_partial[np.logical_not(features_partial.duplicated())]
features_partial

# %%
# save feature full
features_partial.to_csv(features_partial_path, sep='\t', index=False)

# %%
# rename columns
features_full_dic = {
    "feature": features, 
    "intogen": features_intogen,
    "absplice": features_absplice,
    "outlier": features_outlier,
    "coess_cluster": features_coess_modules,
    "joint_embedding_1": features_joint_embedding_1,
    }        

# %%
features_full_list = list(features_full_dic.values())
features_full = reduce(lambda  left,right: pd.merge(left,right,on='gene_id', how='left'), features_full_list).fillna(0)
features_full = features_full.loc[:, (features_full != 0).any(axis=0)]
features_full = features_full[np.logical_not(features_full.duplicated())]
features_full

# %%
# save feature full
features_full.to_csv(features_full_path, sep='\t', index=False)

# %%
