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

# %% [markdown]
# # filter criteria
# - pr coding gene
# - cutoff (can be defined in code chunk 6)
# - tissue by cohort_group
# - merge by gene

# %%
# import modules
import os
import sys
import re
import glob
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm

# %%
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/figure_3_prep_absplice.p"
pickle.dump(snakemake, open(snakemake_path, "wb")) 


# %%
# file = open("/s/project/vale/driver_prediction_202304/processed_data/snakemake/figure_3_prep_absplice.p",'rb')
# snakemake = pickle.load(file)

# %%
def add_cancer_driver_info(res, cgc_cancer_gene, cgc_leukemia_gene, cgc_AML_gene, intogen_cancer_gene):
    
    res['isCGC'] = np.isin(res['gene_id'], cgc_cancer_gene['ENSGid'])
    res['isLeukemia'] = np.isin(res['gene_id'], cgc_leukemia_gene['ENSGid'])
    res['isAML'] = np.isin(res['gene_id'], cgc_AML_gene['ENSGid'])
    res['isOCG'] = np.isin(res['gene_id'], cgc_cancer_gene.ENSGid[cgc_cancer_gene.RoleinCancer.str.contains('oncogene').fillna(False)])
    res['isTSG'] = np.isin(res['gene_id'], cgc_cancer_gene.ENSGid[cgc_cancer_gene.RoleinCancer.str.contains('TSG').fillna(False)])
    
    res = res.merge(cgc_cancer_gene[['ENSGid', 'RoleinCancer']], left_on='gene_id', right_on='ENSGid', how='left').rename(columns={'RoleinCancer':'RoleCGC'})

    intogen_cancer_gene = intogen_cancer_gene.rename(columns={'SYMBOL':'gene_name'})
    res['isIntOGen'] = np.isin(res['gene_name'], intogen_cancer_gene['gene_name'])
    res = res.merge(intogen_cancer_gene[['gene_name', 'ROLE']].drop_duplicates(), on='gene_name', how='left').rename(columns={'ROLE':'RoleIntogen'})    
    
    return res


# %%
cgc_cancer_gene_path = snakemake.params.cgc_cancer_gene_processed
intogen_cancer_gene_path = snakemake.params.intogen_cancer_gene
absplice_res = snakemake.input.absplice_res
absplice_res_var = snakemake.input.absplice_res_var
out_path = snakemake.output.absplice
out_var_path = snakemake.output.absplice_var

# %%
cgc_cancer_gene = pd.read_csv(cgc_cancer_gene_path, sep='\t')
cgc_leukemia_gene = cgc_cancer_gene[["L" in x for x in cgc_cancer_gene['TissueType']]]
cgc_AML_gene = cgc_cancer_gene[["AML" in str(x) for x in cgc_cancer_gene['TumourTypes(Somatic)']]]
# cgc_cancer_gene.head()

intogen_cancer_gene = pd.read_csv(intogen_cancer_gene_path, sep='\t')
# intogen_cancer_gene.head()

# %% [markdown]
# # gene level

# %%
absplice_gene_df = pd.read_csv(absplice_res, sep='\t')

# %%
absplice_gene_df = add_cancer_driver_info(absplice_gene_df, cgc_cancer_gene, cgc_leukemia_gene, cgc_AML_gene, intogen_cancer_gene)
absplice_gene_df = absplice_gene_df.drop('ENSGid', axis=1)

# %%
absplice_gene_df

# %%
absplice_gene_df.to_csv(out_path, sep='\t', index=False)

# %% [markdown]
# # var level

# %%
absplice_var_df = pd.read_csv(absplice_res_var, sep='\t')

# %%
absplice_var_df = add_cancer_driver_info(absplice_var_df, cgc_cancer_gene, cgc_leukemia_gene, cgc_AML_gene, intogen_cancer_gene)
absplice_var_df = absplice_var_df.drop('ENSGid', axis=1)

# %%
absplice_var_df

# %%
absplice_var_df.to_csv(out_var_path, sep='\t', index=False)
