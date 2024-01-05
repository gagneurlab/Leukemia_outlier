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
import pandas as pd
import numpy as np
import pickle

# %%
from scipy.stats import fisher_exact
import statsmodels.stats.multitest as smt

# %%
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/figure_4_prep_enrichment.p"
pickle.dump(snakemake, open(snakemake_path, "wb")) 

# %%
# # read in snakemake object
# file = open("/s/project/vale/driver_prediction_202304/processed_data/snakemake/figure_4_prep_enrichment.p",'rb')
# snakemake = pickle.load(file)

# %%
cgc_cancer_gene = pd.read_csv(snakemake.params.cgc_cancer_gene_processed, sep='\t')
cgc_leukemia_gene = cgc_cancer_gene[["L" in x for x in cgc_cancer_gene['TissueType']]]
cgc_AML_gene = cgc_cancer_gene[["AML" in str(x) for x in cgc_cancer_gene['TumourTypes(Somatic)']]]
# cgc_cancer_gene.head()

samp_anno = pd.read_csv(snakemake.params.samp_anno, sep='\t')
samp_anno = samp_anno.rename(columns={'ArrayID':'sampleID'})
samp_anno_explode = samp_anno.assign(DROP_GROUP=samp_anno.DROP_GROUP.str.split(',')).explode('DROP_GROUP')

gencode = pd.read_csv(snakemake.params.gencode)
gencode = gencode.rename(columns={'gene_id':'gene_id_unique', 'gene_name':'geneSymbol'})
gencode['gene_id'] = [i.split('.', 1)[0] for i in gencode.gene_id_unique]

# %% [markdown]
# # define input list

# %%
category_to_check = 'Diag'

sample_group = snakemake.params.single_group[0]
sample_group_id = samp_anno_explode.loc[samp_anno_explode.DROP_GROUP==sample_group, 'sampleID'].unique()

# read in sa, subset for activation analysis samples
samp_anno = samp_anno[[sample_group in str(x) for x in samp_anno['DROP_GROUP']]]
samp_anno['category'] = samp_anno[[category_to_check]]
samp_anno_enrich = samp_anno[['sampleID', 'category']]

# samp_anno_enrich shoud contain 2 columns: sampleID and category
samp_anno_enrich

# %%
input_res_list = ["ac", "or_up", "or_dn", "fr", "absplice"]

# %%
category_fisher_all_df = pd.DataFrame(columns=['gene_id', 'Diag', 'oddsr', 'p_greater', 'p_greater_adjust',
       'geneSymbol', 'isCGC', 'RoleinCancer', 'TumourTypes(Somatic)', 'input_res'])

for input_res in input_res_list:

    if input_res == "ac":
        # read in activation result
        ac_res = pd.read_csv(snakemake.input.ac_res, sep='\t')
        ac_res = ac_res.loc[ac_res.sampleID.isin(sample_group_id)]

        ac_res = ac_res.rename(columns={'geneID':'gene_id_unique'})
        ac_res['gene_id'] = [i.split('.', 1)[0] for i in ac_res.gene_id_unique]
        ac_res = ac_res.rename(columns={'hgncSymbol':'geneSymbol'})

        ac_res = ac_res.loc[ac_res.AberrantByGene > 1]

        res_en = ac_res[['gene_id', 'geneSymbol', 'sampleID', 'aberrant']]
        res_en = res_en.fillna(value={'aberrant':False})

    if input_res == "or_up":
        # read in OUTRIDER result
        or_res = pd.read_csv(snakemake.input.or_res, sep='\t')
        or_res = or_res.loc[or_res.sampleID.isin(sample_group_id)]

        or_res = or_res.rename(columns={'geneID':'gene_id_unique'})
        or_res['gene_id'] = [i.split('.', 1)[0] for i in or_res.gene_id_unique]
        or_res = or_res.rename(columns={'hgncSymbol':'geneSymbol'})

        or_res['AberrantByGene'] = or_res['gene_id'].map(or_res['gene_id'].value_counts())
        or_res = or_res.loc[or_res.zScore > 0]
        or_res = or_res.loc[or_res.AberrantByGene > 1]

        res_en = or_res[['gene_id', 'geneSymbol', 'sampleID', 'aberrant']]
        res_en = res_en.fillna(value={'aberrant':False})

    if input_res == "or_dn":
        # read in OUTRIDER result
        or_res = pd.read_csv(snakemake.input.or_res, sep='\t')
        or_res = or_res.loc[or_res.sampleID.isin(sample_group_id)]

        or_res = or_res.rename(columns={'geneID':'gene_id_unique'})
        or_res['gene_id'] = [i.split('.', 1)[0] for i in or_res.gene_id_unique]
        or_res = or_res.rename(columns={'hgncSymbol':'geneSymbol'})

        or_res['AberrantByGene'] = or_res['gene_id'].map(or_res['gene_id'].value_counts())
        or_res = or_res.loc[or_res.zScore < 0]
        or_res = or_res.loc[or_res.AberrantByGene > 1]

        res_en = or_res[['gene_id', 'geneSymbol', 'sampleID', 'aberrant']]
        res_en = res_en.fillna(value={'aberrant':False})

    if input_res == "fr":
        # read in FRASER result
        fr_res = pd.read_csv(snakemake.input.fr_res, sep='\t')
        fr_res = fr_res.loc[fr_res.sampleID.isin(sample_group_id)]
        fr_res = fr_res.rename(columns={'hgncSymbol':'geneSymbol'})
        fr_res = pd.merge(fr_res, gencode[['gene_id', 'geneSymbol']], on='geneSymbol', how='left')

        fr_res = pd.merge(fr_res, fr_res.groupby('gene_id', as_index=False).apply(lambda x: len(x))).rename(columns={None:'AberrantByGene'})
        fr_res = fr_res.loc[fr_res.AberrantByGene > 1]
        fr_res = fr_res.assign(aberrant = True)

        res_en = fr_res[['gene_id', 'geneSymbol', 'sampleID', 'aberrant']]

    if input_res == "absplice":
        # read in abSplice result
        absplice_res = pd.read_csv(snakemake.input.absplice_res, sep='\t')
        absplice_res = absplice_res.loc[absplice_res.sampleID.isin(sample_group_id)]
        absplice_res.loc[absplice_res['AbSplice_DNA'] >= 0.2, 'aberrant'] = True
        absplice_res.loc[absplice_res['AbSplice_DNA'] < 0.2, 'aberrant'] = False
        absplice_res = absplice_res.rename(columns={"gene_name":"geneSymbol"})

        absplice_res = pd.merge(absplice_res, absplice_res.groupby('gene_id', as_index=False).apply(lambda x: len(x))).rename(columns={None:'AberrantByGene'})
        absplice_res = absplice_res.loc[absplice_res.AberrantByGene > 1]

        res_en = absplice_res[['gene_id', 'geneSymbol', 'sampleID', 'aberrant']]

        # res_en should contain 4 columns: gene_id, geneSymbol, sampleID, aberrant
        # res_en

    # fisher test
    # initialize result table
    category_fisher_df = pd.DataFrame(columns=['gene_id','Diag','oddsr','p_greater'])

    for gene_to_check in res_en.gene_id.unique():

        res_en_gene = res_en.loc[res_en.gene_id==gene_to_check]
        res_en_gene_all = pd.merge(res_en_gene, samp_anno_enrich[['sampleID', 'category']], on='sampleID', how='outer')
        res_en_gene_all = res_en_gene_all.fillna(value={'aberrant':False})

        diag_list = res_en_gene_all.loc[res_en_gene_all.aberrant==True].category.unique()

        for diag_to_check in diag_list:
            fisher_df = res_en_gene_all
            fisher_df['inCategory'] = fisher_df.category == diag_to_check
            fisher_df = fisher_df[['aberrant', 'inCategory']]
        #     fisher_df

            con_tab = pd.crosstab(index=fisher_df['aberrant'], columns=fisher_df['inCategory'])
        #     con_tab

            oddsr, p_greater = fisher_exact(con_tab, alternative='greater')

            category_fisher_df_sub = pd.DataFrame(data={'gene_id': [gene_to_check], 'Diag': [diag_to_check], 'oddsr': [oddsr], 'p_greater': [p_greater]}) 
            category_fisher_df = pd.concat([category_fisher_df, category_fisher_df_sub])

    category_fisher_df = category_fisher_df.assign(p_greater_adjust = 
                                                   smt.multipletests(category_fisher_df.p_greater, method='fdr_bh')[1])
    category_fisher_df = category_fisher_df.sort_values(by=['p_greater_adjust', 'oddsr'], ascending=True)

    category_fisher_df = pd.merge(category_fisher_df, res_en[['gene_id', 'geneSymbol']].drop_duplicates(), on='gene_id', how='left')
    category_fisher_df['isCGC'] = np.isin(category_fisher_df.gene_id, cgc_cancer_gene.ENSGid)
    category_fisher_df = pd.merge(category_fisher_df, cgc_cancer_gene[['ENSGid', 'RoleinCancer', 'TumourTypes(Somatic)']], left_on='gene_id', right_on='ENSGid', how='left')
    category_fisher_df = category_fisher_df.drop('ENSGid', axis=1)
    category_fisher_df = category_fisher_df.assign(input_res=input_res)
    
    print(input_res + " done")

    category_fisher_all_df = pd.concat([category_fisher_all_df, category_fisher_df])

# %%
category_fisher_all_df

# %%
category_fisher_all_df.to_csv(snakemake.output.diag_fisher, sep='\t', index=False)
