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

# %% [markdown]
# Author: Kim Anh Lilian Le, Xueqi Cao

# %%
import os
import pickle
import pandas as pd
import numpy as np

# %%
snakemake_path = snakemake.params.projectPath + "/processed_data/snakemake/generate_coess_feature.p"
pickle.dump(snakemake, open(snakemake_path, "wb"))

# %%
# file = open("/s/project/vale/driver_prediction_202207/processed_data/snakemake/generate_coess_feature.p",'rb')
# snakemake = pickle.load(file)

# %% [markdown]
# # Process gene code 

# %%
coess_modules_path = snakemake.input.coess_modules
joint_embedding_1_path = snakemake.input.joint_embedding_1
gencode_path = snakemake.params.gencode
out_dir = os.path.dirname(snakemake.output.coess_modules)
joint_embedding_1_out_path = snakemake.output.joint_embedding_1

# %%
# coess_modules_path = "/s/project/tcga/driver_prediction_model/coess_cluster/coess_supplementary_data/coessential_modules.csv"
# gencode_path = "/s/project/mll/drop_2020apr/processed_data/aberrant_expression/v33b/gene_name_mapping_v33b.tsv"
# out_path = '/s/project/vale/driver_prediction_202204/processed_data/coess_feature'

# %%
# create output dir if not exist
if not os.path.isdir(out_dir): 
    os.mkdir(out_dir)

# %%
gencode = pd.read_csv(gencode_path)

gencode = gencode.rename(columns={'gene_id':'gene_id_unique'})
gencode['gene_id'] = [i.split('.', 1)[0] for i in gencode.gene_id_unique]

features = gencode[["gene_name", "gene_id"]]

# %% [markdown]
# # Co-essential modules

# %%
coess_sup_df = pd.read_csv(coess_modules_path, sep = '\t', skiprows=2, low_memory=False)
# coess_sup_df

# %%
d_range = coess_sup_df["d"].unique()
d_range

# %%
for d in d_range:
    sub_df = coess_sup_df[coess_sup_df["d"] == d]
    outPath = out_dir + f"/coess_cluster-d{d}.tsv"

    # Split for only the row where genes are
    index_split = sub_df.columns.tolist().index("Genes")
    last_index = len(sub_df.columns.tolist())

    filt_df = sub_df.iloc[:, index_split:last_index]
    
    # Get list of genes, remove NAs
    out_ls = []

    for l in filt_df.values:
        cleaned_list = [x for x in l if not pd.isnull(x)]
        out_ls.append(cleaned_list)
                
    out_df = pd.DataFrame(columns=['gene_name', f'coess_cluster-d{d}'])

    for i, gene_ls in enumerate(out_ls, 1):
        cluster_ls = [i] * len(gene_ls)
        sub_df = pd.DataFrame({'gene_name': gene_ls, f'coess_cluster-d{d}': cluster_ls})
        out_df = pd.concat([out_df, sub_df])
        
    out_df = pd.merge(features, out_df, on=["gene_name"], how='right')
    
    # dcast the df
    out_df = out_df.pivot_table(index=['gene_id', 'gene_name'], columns=[f"coess_cluster-d{d}"], aggfunc=len).fillna(0).reset_index()
    out_df = out_df.drop('gene_name', axis=1)
    out_df.columns = out_df.columns[:1].tolist() + [f'coess_cluster-d{d}-' + str(col) for col in out_df.columns[1:]]
    
    # Save df
    print(f"Saving in {outPath}")
    out_df.to_csv(outPath, sep='\t', index=False)
    
    # Nr. unique genes per module
    print(f"d = {d}")
    print(f"Nr. unique genes: {len(out_df)}")
    print(f"Nr. cluster in: {len(out_df.columns)-1} \n")

# %% [markdown]
# # joint embedding 1

# %%
joint_embedding_1 = pd.read_csv(joint_embedding_1_path, sep="\t")
joint_embedding_1.columns = joint_embedding_1.columns[:1].tolist() + [f'joint_embedding_1-' + str(col) for col in joint_embedding_1.columns[1:]]

print(f"Saving in {joint_embedding_1_out_path}")
joint_embedding_1.to_csv(joint_embedding_1_out_path, sep="\t", index=False)

print(f"Nr. unique genes: {len(joint_embedding_1)}")
print(f"Nr. dimensions in: {len(joint_embedding_1.columns)-1} \n") 
