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
import gseapy as gp
from tqdm import tqdm
from gseapy.plot import gseaplot


# %% [markdown]
# | {es: enrichment score,
# |  NES: normalized enrichment score,
# |  p: P-value,
# |  FDRq_val: FDRq_val,
# |  size: gene set size,
# |  matched_size: genes matched to the data,
# |  genes: gene names from the data set
# |  ledge_genes: leading edge genes}

# %%
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/figure_5_GSEA.p"
pickle.dump(snakemake, open(snakemake_path, "wb")) 

# %%
# # read in snakemake object
# file = open("/s/project/vale/driver_prediction_202303_2/processed_data/snakemake/figure_5_GSEA.p",'rb')
# snakemake = pickle.load(file)

# %%
gencode = pd.read_csv(snakemake.params.gencode)
gencode = gencode.loc[gencode.gene_type=='protein_coding']

# %%
output_path = snakemake.output.diag_fisher_gsea
output_dir = snakemake.params.projectPath + "/manuscript/figure_5/plot_data/gsea"
if not os.path.exists(output_dir):
   os.makedirs(output_dir)

# %%
hallmark_gs = gp.get_library(name="MSigDB_Hallmark_2020")
# print(hallmark_gs['Myc Targets V1'])
# print(len(hallmark_gs['Myc Targets V1']))

# %% [markdown]
# # GSEA

# %%
diag_fisher = pd.read_csv(snakemake.input.diag_fisher, sep='\t')
# diag_fisher = pd.read_csv('/s/project/vale/driver_prediction_202303_2/manuscript/figure_5/plot_data/diag_fisher.tsv', sep='\t')

# %%
gene_sets_list = ['MSigDB_Hallmark_2020', 'MSigDB_Oncogenic_Signatures']

# %%
test_diag_list = diag_fisher.Diag.value_counts()[diag_fisher.Diag.value_counts()>1].index
test_diag_list

# %%
tool_list = diag_fisher.input_res.unique()
tool_list

# %%
product_list = list(itertools.product(tool_list, test_diag_list, gene_sets_list))
len(product_list)

# %%
for i in tqdm(product_list):
    tool_selected = i[0]
    test_diag = i[1]
    gene_sets_selectecd = i[2]

    diag_fisher_sub = diag_fisher.loc[diag_fisher.input_res==tool_selected]
    diag_fisher_sub = diag_fisher_sub.loc[diag_fisher_sub.Diag==test_diag]
    # diag_fisher_sub = diag_fisher_sub.query('p_greater_adjust<0.05')
    diag_fisher_sub
    
    print(len(diag_fisher_sub. index))
    if len(diag_fisher_sub. index)>=2: 

        gsea_input = diag_fisher_sub[['geneSymbol', 'oddsr', 'p_greater_adjust']]
        # saturate Inf with max
        # gsea_input.loc[gsea_input.oddsr==np.Inf, 'oddsr'] = max(gsea_input.loc[gsea_input.oddsr!=np.Inf]['oddsr'])
        gsea_input = gsea_input.merge(gencode[['gene_name']].rename(columns={'gene_name':'geneSymbol'}), on='geneSymbol', how='right')
        gsea_input = gsea_input.fillna({'oddsr':0,'p_greater_adjust':1})

        # use oddsr
        gsea_input = gsea_input.rename(columns={'geneSymbol':0, 'oddsr':1}).sort_values(1, ascending=False).drop(['p_greater_adjust'], axis=1)

        # use p_greater_adjust
        # gsea_input = gsea_input.rename(columns={'geneSymbol':0, 'p_greater_adjust':1}).sort_values(1, ascending=False).drop(['oddsr'], axis=1)

        # Run GSEA pre ranked with hallmarks_genes
        gsea_obj = gp.prerank(rnk=gsea_input, gene_sets=gene_sets_selectecd,
                             processes=10,
                             permutation_num=1000, # reduce number to speed up testing
                             outdir=output_dir, format='png', seed=5)

        # save gsea
        gsea_obj_path = os.path.join(output_dir, 'gsea_obj' + '_' + gene_sets_selectecd + '_' + test_diag.replace("/", "_") + ".p")
        pickle.dump(gsea_obj, open(gsea_obj_path, "wb"))

        gsea_res = gsea_obj.res2d
        gsea_res.columns = [i.replace(" ", "") for i in gsea_res.columns]
        gsea_res.columns = [i.replace("-", "_") for i in gsea_res.columns]
        gsea_res = gsea_res.assign(NES_abs = gsea_res.NES.abs())
        gsea_res = gsea_res.assign(input_res=tool_selected)
        gsea_res = gsea_res.assign(test_diag=test_diag)
        gsea_res = gsea_res.assign(gsea_list=gene_sets_selectecd)
        gsea_res = gsea_res.assign(obj_path=gsea_obj_path)

        if 'gsea_res_all' not in locals():
            gsea_res_all = gsea_res
        else:
            gsea_res_all = pd.concat([gsea_res_all, gsea_res])

# %%
gsea_res_all

# %%
gsea_res_all.to_csv(output_path, sep='\t', index=False)

# %% [markdown]
# # visualization

# %%
# gsea_res_all = pd.read_csv(output_path, sep='\t')
gsea_res_all

# %%
# Check table with the enrichment scores
gsea_res_sig = gsea_res_all.sort_values(by="NES_abs", ascending=False).query('FDRq_val<0.1').reset_index()
gsea_res_sig

# %%
gsea_res_sig.Term.value_counts()

# %%
viz_num = 0

# %%
# Run gsea plot
term_viz = gsea_res_sig.Term[viz_num]
file = open(gsea_res_sig.obj_path[viz_num],'rb')
gsea_obj_viz = pickle.load(file)

gseaplot(rank_metric=gsea_obj_viz.ranking, term=term_viz, 
         **gsea_obj_viz.results[term_viz])

# %%
# Check table with the enrichment scores
gsea_res_sig.loc[gsea_res_sig.input_res=='ac']

    # %%
    tool_selected = 'ac'
    test_diag = 'MDS/MPN-RS-T'

    diag_fisher_sub = diag_fisher.loc[diag_fisher.input_res==tool_selected]
    diag_fisher_sub = diag_fisher_sub.loc[diag_fisher_sub.Diag==test_diag]
    # diag_fisher_sub = diag_fisher_sub.query('p_greater_adjust<0.05')
    diag_fisher_sub

