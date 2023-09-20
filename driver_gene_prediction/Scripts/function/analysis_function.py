import os
import re
import pandas as pd
import numpy as np
import math


def simplify_exp(exp_viz):
    # exp_viz = exp_viz.assign(label_gene_list = exp_viz.apply(lambda x: re.sub('CGC_', '', x['label_gene_list']), axis=1))
    exp_viz.loc[exp_viz.sample_group == 'leukemia_14group' , 'sample_group'] = 'leu14'
    exp_viz.loc[exp_viz.intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' , 'intogen_input_feature'] = 'g7tools'
    exp_viz.loc[exp_viz.coess_input_feature == 'coess_cluster' , 'coess_input_feature'] = 'coessC'
    
    return exp_viz


def add_exp_tag(exp_viz, exp_tag_common_diaplay, exp_tag_individual_display_list):

    exp_viz_sub = exp_viz[['sample_group', 'label_gene_list', 'model_method',
                           'intogen_input_feature', 'outlier_input_feature', 'coess_input_feature', ]]

    for col in exp_viz_sub:
        if(len(exp_viz_sub.loc[:,col].unique()) != 1):
            exp_viz_sub.pop(col)
            
    exp_tag_common = "-".join(exp_viz_sub.iloc[0])
    exp_viz = exp_viz.assign(exp_tag_common = exp_tag_common)
    
    exp_viz_sub = exp_viz[['sample_group', 'label_gene_list', 'model_method',
                           'intogen_input_feature', 'outlier_input_feature', 'coess_input_feature', ]]

    for col in exp_viz_sub:
        if(len(exp_viz_sub.loc[:,col].unique()) == 1):
            exp_viz_sub.pop(col)
            
    exp_viz = exp_viz.assign(exp_tag_individual = exp_viz_sub.apply(lambda x: "-".join(x.replace("", np.nan).dropna()), axis=1))

    if exp_tag_common_diaplay != False:
        exp_viz['exp_tag_common_diaplay'] = exp_tag_common_diaplay
    else:
        exp_tag_common_diaplay = ": " + exp_tag_common
        exp_viz['exp_tag_common_diaplay'] = exp_tag_common_diaplay
        
    if exp_tag_individual_display_list != False:
        exp_viz['exp_tag_individual_display'] = exp_tag_individual_display_list
    else:
        exp_viz['exp_tag_individual_display'] = exp_viz['exp_tag_individual']
        exp_tag_individual_display_list = exp_viz['exp_tag_individual_display'].to_list()

    return exp_viz, exp_tag_common, exp_tag_common_diaplay, exp_tag_individual_display_list


def add_result_paths(exp_viz, project_dir, intogen_dir):
    exp_viz = exp_viz.assign(res_pre_path = exp_viz.apply(lambda x: project_dir + '/processed_results/pre/result_' + str(x['experiment_no']) + '.tsv', axis=1))
    exp_viz = exp_viz.assign(coeff_path = exp_viz.apply(lambda x: project_dir + '/processed_results/pre/coeff_' + str(x['experiment_no']) + '.tsv', axis=1))
    exp_viz = exp_viz.assign(res_post_path = exp_viz.apply(lambda x: project_dir + '/processed_results/post/result_' + str(x['experiment_no']) + '.tsv', axis=1))
    exp_viz = exp_viz.assign(rp_path = exp_viz.apply(lambda x: project_dir + '/processed_results/post/rp_' + str(x['experiment_no']) + '.tsv', axis=1))
    exp_viz = exp_viz.assign(vep_path = exp_viz.apply(lambda x: intogen_dir + '/steps/vep/MLL_WGS_MLL_' + str(x['sample_group']).upper() + '.tsv.gz', axis=1))
    exp_viz = exp_viz.assign(dndscv_path = exp_viz.apply(lambda x: intogen_dir + '/steps/dndscv/MLL_WGS_MLL_' + str(x['sample_group']).upper() + '.dndscv.tsv.gz', axis=1))

    return exp_viz

def set_role(data, distance_threshold=0.1):
    """Set the role according to the DNDS output"""
    if data['wmis_cv'] < 1 and data['wnon_cv'] < 1:  # threshold
        return "ambiguous"
    
    # Check wmis
    wmis = data['wmis_cv']
    if wmis >= 1 and data["n_mis"] == 0:
        wmis = 1

    # Check wnon
    wnon = data['wnon_cv']
    if wnon >= 1 and data["n_non"] == 0:
        wnon = 1
        
    # Those cases with w_non and w_mis <=1 are not informative
    if wnon <= 1 and wmis <= 1:
        return "ambiguous"

    distance = (wmis - wnon) / math.sqrt(2)
    if distance_threshold is not None and abs(distance) < distance_threshold:
        return "ambiguous"
    else:
        if distance > 0:
            return 'Act'
        elif distance < 0:
            return 'LoF'
        else:
            return "ambiguous"
        
def estimate_role(data):
    distance = data['CGC_leukemia_TSG'] - data['CGC_leukemia_OCG']

    if distance > 10:
        return 'Act'
    elif distance < -10:
        return 'LoF'
    else:
        return "ambiguous"