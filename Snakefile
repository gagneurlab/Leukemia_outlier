import pathlib
import wbuild
import glob
import pandas as pd
import numpy as np
from pathlib import Path

config['wBuildPath'] =  str(pathlib.Path(wbuild.__file__).parent)

configfile: "wbuild.yaml"
include: config['wBuildPath'] + "/wBuild.snakefile"

htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"


if not os.path.isdir(config['tmpDirPath']): 
  os.mkdir(config['tmpDirPath'])

if not os.path.isdir(config['projectPath'] + "/processed_data"): 
  os.mkdir(config['projectPath'] + "/processed_data")
  os.mkdir(config['projectPath'] + "/processed_data/snakemake")
  
if not os.path.isdir(config['projectPath'] + "/processed_results"): 
  os.mkdir(config['projectPath'] + "/processed_results")
  os.mkdir(config['projectPath'] + "/processed_results/pre")
  os.mkdir(config['projectPath'] + "/processed_results/post")
  
if not os.path.isdir(config['projectPath'] + "/Output"): 
  os.mkdir(config['projectPath'] + "/Output")
  os.mkdir(config['projectPath'] + "/Output/html")
  os.mkdir(config['projectPath'] + "/Output/plot")

        
experimentDesign = pd.read_csv(config["experimentDesign"], sep='\t')
experimentDesign = experimentDesign.fillna('')
#experimentDesign = experimentDesign.loc[
#    (np.isin(experimentDesign.model_method, ['rf']))
#]


sampleAnnotation = pd.read_csv(config["sampleAnnotation"], sep='\t')
absplice_files = sampleAnnotation.DNA_ABSPLICE_FILE.dropna()

rule all:
    input: 
        rules.Index.output, 
        htmlOutputPath + "/readme.html",
        config['tmpDirPath'] + "/generate_label.done",
        config['tmpDirPath'] + "/prediction.done",
        config['tmpDirPath'] + "/postprocessing.done"
    output: 
        touch(config['tmpDirPath'] + "/all.done")
        
rule extract_intogen_score:
    input:
        inputDir=config['intogenDir']
    output:
        score=config["projectPath"] + '/processed_data/intogen_feature/score/cbase.tsv'
    params:
        projectPath=config['projectPath']
    script: 
        "Scripts/preprocessing/intogen_feature/extract_intogen_score.py"
       
rule generate_intogen_feature:
    input:
        score=config["projectPath"] + '/processed_data/intogen_feature/score/cbase.tsv'
    output:
        feature=config["projectPath"] + '/processed_data/intogen_feature/feature/cbase.tsv'
    params:
        projectPath=config['projectPath'],
        intogenDir=config['intogenDir'],
        gencode=config['gencode'],
        single_group=config['cohort']['single_group'],
        samp_anno=config['sampleAnnotation']
    script: 
        "Scripts/preprocessing/intogen_feature/generate_intogen_feature.py"
        
rule generate_outlier_feature:
    input:
        rules.preprocessing_outlier_feature_generate_outrider_features_R.output,
        rules.preprocessing_outlier_feature_generate_outrider_features_groupwise_R.output,
        rules.preprocessing_outlier_feature_generate_activation_features_R.output,
        rules.preprocessing_outlier_feature_generate_activation_features_groupwise_R.output,
        rules.preprocessing_outlier_feature_generate_fraser_features_groupwise_R.output
    output:
        touch(config['tmpDirPath'] + "/generate_outlier_feature.done")


rule read_in_absplice_res:
    input:
        absplice=absplice_files
    output:
        absplice_res=config["projectPath"] + '/processed_data/absplice_res.tsv',
        absplice_res_var=config["projectPath"] + '/processed_data/absplice_res_var.tsv'
    params:
        single_group=config['cohort']['single_group'],
        projectPath=config['projectPath'],
        gencode=config['gencode'],
        samp_anno=config['sampleAnnotation'],
    conda:
        "envs/absplice.yaml"
    script: 
        "Scripts/preprocessing/absplice_feature/read_in_absplice_res.py"
        
        
rule generate_absplice_feature:
    input:
        absplice_res=config["projectPath"] + '/processed_data/absplice_res.tsv'
    output:
        absplice=config["projectPath"] + '/processed_data/absplice_feature/absplice.tsv',
    params:
        samp_groups=config['cohort']['groups'],
        single_group=config['cohort']['single_group'],
        projectPath=config['projectPath'],
        gencode=config['gencode'],
        samp_anno=config['sampleAnnotation'],
    script: 
        "Scripts/preprocessing/absplice_feature/generate_absplice_feature.py"
        
        
rule generate_coess_feature:
    input:
        coess_modules=config['coess_modules'],
        joint_embedding_1=config['joint_embedding_1']
    output:
        coess_modules=config["projectPath"] + '/processed_data/coess_feature/coess_cluster-d0.9.tsv',
        joint_embedding_1=config["projectPath"] + '/processed_data/coess_feature/joint_embedding_1.tsv'
    params:
        projectPath=config['projectPath'],
        gencode=config['gencode']
    script: 
        "Scripts/preprocessing/coess_feature/generate_coess_feature.py"
        
rule merge_feature:
    resources:
        mem_mb = 8000
    input:
        config['tmpDirPath'] + "/generate_outlier_feature.done",
        intogen=config["projectPath"] + '/processed_data/intogen_feature/feature/cbase.tsv',
        absplice=config["projectPath"] + '/processed_data/absplice_feature/absplice.tsv',
        coess_modules=config["projectPath"] + '/processed_data/coess_feature/coess_cluster-d0.9.tsv',
        joint_embedding_1=config["projectPath"] + '/processed_data/coess_feature/joint_embedding_1.tsv'
    output:
        features_partial=config["projectPath"] + '/processed_data/features_partial.tsv',
        features_full=config["projectPath"] + '/processed_data/features_full.tsv'
    params:
        outlierDir=config["projectPath"] + '/processed_data/outlier_feature',
        projectPath=config['projectPath'],
        gencode=config['gencode']
    script: 
        "Scripts/preprocessing/merge_feature.py"
        
rule generate_label:
    input:
        rules.preprocessing_label_generate_label_list_R.output
    output:
        touch(config['tmpDirPath'] + "/generate_label.done")

rule prediction:
    resources:
        mem_mb = 6000
    params:
        sampID = "{exp_no}",
        projectPath=config['projectPath']
    input:
        label=rules.preprocessing_label_generate_label_list_R.output,
        experimentDesign=config["projectPath"] + '/experiment_design.tsv',
        features_partial=config["projectPath"] + '/processed_data/features_partial.tsv',
        features_full=config["projectPath"] + '/processed_data/features_full.tsv'
    output:
        result=config["projectPath"] + '/processed_results/pre/result_{exp_no}.tsv',
        coeff=config["projectPath"] + '/processed_results/pre/coeff_{exp_no}.tsv'
    script: 
        "Scripts/prediction/prediction.py"
        
rule prediction_done:
    input:
        result_post=expand(config["projectPath"] + '/processed_results/pre/result_{exp_no}.tsv', exp_no=experimentDesign['experiment_no']),
        coeff=expand(config["projectPath"] + '/processed_results/pre/coeff_{exp_no}.tsv', exp_no=experimentDesign['experiment_no'])
    output:
        touch(config['tmpDirPath'] + "/prediction.done")
        
rule postprocessing:
    params:
        sampID = "{exp_no}",
        tmpDirPath=config["tmpDirPath"],
        projectPath = config["projectPath"],
        intogenDir = config["intogenDir"],
        outriderDir = config["outriderDir"],
        fraserDir = config["fraserDir"],
        absplice = config["abspliceDir"],
        coess_modules = config["coess_modules"],
        samp_anno=config['sampleAnnotation'],
        gencode=config['gencode'],
        intogen_cancer_gene=config["IntOGen_cancer_gene"],
        cgc_cancer_gene_processed=rules.preprocessing_label_generate_label_list_R.output[0],
        predictedConsequence = config["predictedConsequence"]
    input:
        config['tmpDirPath'] + "/prediction.done",
        experimentDesign=config["projectPath"] + '/experiment_design.tsv',
        features_partial=config["projectPath"] + '/processed_data/features_partial.tsv',
        features_full=config["projectPath"] + '/processed_data/features_full.tsv',
        result=config["projectPath"] + '/processed_results/pre/result_{exp_no}.tsv'
    output:
        result_post=config["projectPath"] + '/processed_results/post/result_{exp_no}.tsv',
        rank_proportion=config["projectPath"] + '/processed_results/post/rp_{exp_no}.tsv'
    script: 
        "Scripts/postprocessing/postprocessing.py"
        
rule postprocessing_done:
    input:
        result_post=expand(config["projectPath"] + '/processed_results/post/result_{exp_no}.tsv', exp_no=experimentDesign['experiment_no']),
        rank_proportion=expand(config["projectPath"] + '/processed_results/post/rp_{exp_no}.tsv', exp_no=experimentDesign['experiment_no'])
    output:
        touch(config['tmpDirPath'] + "/postprocessing.done")
        
rule absplice_tab:
    input:
        absplice_res=config["projectPath"] + '/processed_data/absplice_res.tsv'
    output:
        absplice_tab=config["projectPath"] + '/manuscript/sup_table/absplice_tab.csv'
    params:
        single_group=config['cohort']['single_group'],
        projectPath=config['projectPath'],
        gencode=config['gencode'],
        samp_anno=config['sampleAnnotation'],
    script: 
        "Scripts/manuscript/sup_table/absplice_tab.py"
        
rule intogen_tab:
    input:
        inputDir=config['intogenDir']
    output:
        intogen_tab=config["projectPath"] + '/manuscript/sup_table/intogen_tab.csv'
    params:
        projectPath=config['projectPath'],
        manuscriptWording=config["manuscriptWording"]
    script: 
        "Scripts/manuscript/sup_table/intogen_tab.py"
