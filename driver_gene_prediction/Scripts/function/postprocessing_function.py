import os
import pandas as pd
import numpy as np

from prediction_function import *

def aggregate_res_by_gene(res):
    res = pd.pivot_table(res, values='Prediction', index = ['GeneID', 'GeneSymbol', 'Label'], aggfunc=np.median).sort_values(by='Prediction', ascending=False).reset_index()
    res = res.assign(Rank = res.Prediction.rank(ascending=False, method='min'))
    
    return res


def add_cancer_driver_info(res, cgc_cancer_gene, cgc_leukemia_gene, cgc_AML_gene, intogen_cancer_gene):
    res = res.assign(isCGC = np.isin(res['GeneID'], cgc_cancer_gene['ENSGid']))
    res = res.assign(isLeukemia = np.isin(res['GeneID'], cgc_leukemia_gene['ENSGid']))
    res = res.assign(isAML = np.isin(res['GeneID'], cgc_AML_gene['ENSGid']))
    res = res.assign(isOCG = np.isin(res['GeneID'], cgc_cancer_gene.ENSGid[cgc_cancer_gene.RoleinCancer.str.contains('oncogene').fillna(False)]))
    res = res.assign(isTSG = np.isin(res['GeneID'], cgc_cancer_gene.ENSGid[cgc_cancer_gene.RoleinCancer.str.contains('TSG').fillna(False)]))
    res = res.assign(isIntOGen = np.isin(res['GeneSymbol'], intogen_cancer_gene['GeneSymbol']))
    
    res = res.merge(cgc_cancer_gene[['ENSGid', 'RoleinCancer']], left_on='GeneID', right_on='ENSGid', how='left').rename(columns={'RoleinCancer':'RoleCGC'})
    res = res.merge(intogen_cancer_gene[['GeneSymbol', 'ROLE']].drop_duplicates(), on='GeneSymbol', how='left').rename(columns={'ROLE':'RoleIntogen'})
    res = res.drop('ENSGid', axis=1)
    
    return res


def get_intogen_res_no(experiment_no, experiment_df):

    label_gene_list, sample_group, random_seeds, model_method = get_param_training(experiment_no, experiment_df)
    intogen_input_feature, _, _ = get_param_feature(experiment_no, experiment_df)

    experiment_df_sub = experiment_df.loc[(experiment_df.label_gene_list == label_gene_list) *
                      (experiment_df.sample_group == sample_group) *
                      (experiment_df.model_method == model_method)]
    idx = (
        (experiment_df_sub.intogen_input_feature == intogen_input_feature) *
        (experiment_df_sub.outlier_input_feature == "") *
        (experiment_df_sub.coess_input_feature == "")
    )
    
    if any(idx):
        intogen_res_no = experiment_df_sub.loc[idx]['experiment_no'].values[0]
    else:
        intogen_res_no = np.NaN
        
    return(intogen_res_no)


def get_outrider_res_no(experiment_no, experiment_df):
    label_gene_list, sample_group, random_seeds, model_method = get_param_training(experiment_no, experiment_df)
    _, outlier_input_feature, _ = get_param_feature(experiment_no, experiment_df)

    experiment_df_sub = experiment_df.loc[(experiment_df.label_gene_list == label_gene_list) *
                      (experiment_df.sample_group == sample_group) *
                      (experiment_df.model_method == model_method)]

    idx = (
        (experiment_df_sub.intogen_input_feature == "") *
        (experiment_df_sub.outlier_input_feature == outlier_input_feature) *
        (experiment_df_sub.coess_input_feature == "")
    )
    
    if any(idx):
        outrider_res_no = experiment_df_sub.loc[idx]['experiment_no'].values[0]
    else:
        outrider_res_no = np.NaN
    
    return(outrider_res_no)


def get_coess_res_no(experiment_no, experiment_df):
    label_gene_list, sample_group, random_seeds, model_method = get_param_training(experiment_no, experiment_df)
    _, _, coess_input_feature = get_param_feature(experiment_no, experiment_df)

    experiment_df_sub = experiment_df.loc[(experiment_df.label_gene_list == label_gene_list) *
                      (experiment_df.sample_group == sample_group) *
                      (experiment_df.model_method == model_method)]

    idx = (
        (experiment_df_sub.intogen_input_feature == "") *
        (experiment_df_sub.outlier_input_feature == "") *
        (experiment_df_sub.coess_input_feature == coess_input_feature)
    )   
        
    if any(idx):
        coess_res_no = experiment_df_sub.loc[idx]['experiment_no'].values[0]
    else:
        coess_res_no = np.NaN
        
    return(coess_res_no)


def add_intogen_res(res, intogen_input_feature, experiment_no, experiment_df, res_dir):
    if intogen_input_feature != '':
        intogen_res_no = get_intogen_res_no(experiment_no, experiment_df)    

        if np.isnan(intogen_res_no):
            res = res.assign(intogenPred=np.NaN)  
            res = res.assign(intogenPred_Rank=np.NaN) 
        else: 
            intogen_res = pd.read_csv(res_dir + '/result_' + str(intogen_res_no) + '.tsv', sep='\t', index_col=0)
            intogen_res = aggregate_res_by_gene(intogen_res)
            intogen_res = intogen_res.rename(columns={'Prediction':'intogenPred'})
            res = pd.merge(res, intogen_res[['GeneID', 'intogenPred']], on='GeneID')
            res = res.assign(intogenPred_Rank = res.intogenPred.rank(ascending=False, method='min'))
    else:
        res = res.assign(intogenPred=np.NaN)
        res = res.assign(intogenPred_Rank=np.NaN) 
    return res


def add_outrider_res(res, outlier_input_feature, experiment_no, experiment_df, res_dir):
    if outlier_input_feature != '':
        outrider_res_no = get_outrider_res_no(experiment_no, experiment_df)

        if np.isnan(outrider_res_no):
            res = res.assign(outlierPred=np.NaN)     
            res = res.assign(outlierPred_Rank=np.NaN)
        else: 
            outrider_res = pd.read_csv(res_dir + '/result_' + str(outrider_res_no) + '.tsv', sep='\t', index_col=0)
            outrider_res = aggregate_res_by_gene(outrider_res)
            outrider_res = outrider_res.rename(columns={'Prediction':'outlierPred'})
            res = pd.merge(res, outrider_res[['GeneID', 'outlierPred']], on='GeneID')
            res = res.assign(outlierPred_Rank = res.outlierPred.rank(ascending=False, method='min'))
    else:
        res = res.assign(outlierPred=np.NaN)
        res = res.assign(outlierPred_Rank=np.NaN)

    return res


def add_coess_res(res, coess_input_feature, experiment_no, experiment_df, res_dir):
    if coess_input_feature != '':
        coess_res_no = get_coess_res_no(experiment_no, experiment_df)

        if np.isnan(coess_res_no):
            res = res.assign(coessPred=np.NaN)  
            res = res.assign(coessPred_Rank=np.NaN)
        else: 
            coess_res = pd.read_csv(res_dir + '/result_' + str(coess_res_no) + '.tsv', sep='\t', index_col=0)
            coess_res = aggregate_res_by_gene(coess_res)
            coess_res = coess_res.rename(columns={'Prediction':'coessPred'})
            res = pd.merge(res, coess_res[['GeneID', 'coessPred']], on='GeneID')
            res = res.assign(coessPred_Rank = res.coessPred.rank(ascending=False, method='min'))
    else:
        res = res.assign(coessPred=np.NaN)
        res = res.assign(coessPred_Rank=np.NaN)

    return res


def get_wo_intogen_res_no(experiment_no, experiment_df):
    label_gene_list, sample_group, random_seeds, model_method = get_param_training(experiment_no, experiment_df)
    _, outlier_input_feature, coess_input_feature = get_param_feature(experiment_no, experiment_df)

    experiment_df_sub = experiment_df.loc[(experiment_df.label_gene_list == label_gene_list) *
                      (experiment_df.sample_group == sample_group) *
                      (experiment_df.model_method == model_method)]

    idx = (
        (experiment_df_sub.intogen_input_feature == "") *
        (experiment_df_sub.outlier_input_feature == outlier_input_feature) *
        (experiment_df_sub.coess_input_feature == coess_input_feature)
    )
    
    if any(idx):
        intogen_res_no = experiment_df_sub.loc[idx]['experiment_no'].values[0]
    else:
        intogen_res_no = np.NaN

    return(intogen_res_no)


def get_wo_outrider_res_no(experiment_no, experiment_df):
    label_gene_list, sample_group, random_seeds, model_method = get_param_training(experiment_no, experiment_df)
    intogen_input_feature, _, coess_input_feature = get_param_feature(experiment_no, experiment_df)
    
    experiment_df_sub = experiment_df.loc[(experiment_df.label_gene_list == label_gene_list) *
                      (experiment_df.sample_group == sample_group) *
                      (experiment_df.model_method == model_method)]

    idx = (
        (experiment_df_sub.intogen_input_feature == intogen_input_feature) *
        (experiment_df_sub.outlier_input_feature == "") *
        (experiment_df_sub.coess_input_feature == coess_input_feature)
    )
    
    if any(idx):
        outrider_res_no = experiment_df_sub.loc[idx]['experiment_no'].values[0]
    else:
        outrider_res_no = np.NaN
    
    return(outrider_res_no)


def get_wo_coess_res_no(experiment_no, experiment_df):
    label_gene_list, sample_group, random_seeds, model_method = get_param_training(experiment_no, experiment_df)
    intogen_input_feature, outlier_input_feature, _ = get_param_feature(experiment_no, experiment_df)

    experiment_df_sub = experiment_df.loc[(experiment_df.label_gene_list == label_gene_list) *
                      (experiment_df.sample_group == sample_group) *
                      (experiment_df.model_method == model_method)]

    idx = (
        (experiment_df_sub.intogen_input_feature == intogen_input_feature) *
        (experiment_df_sub.outlier_input_feature == outlier_input_feature) *
        (experiment_df_sub.coess_input_feature == "")
    )   
        
    if any(idx):
        coess_res_no = experiment_df_sub.loc[idx]['experiment_no'].values[0]
    else:
        coess_res_no = np.NaN
        
    return(coess_res_no)


def add_wo_intogen_res(res, intogen_input_feature, experiment_no, experiment_df, res_dir):
    if intogen_input_feature != '':
        intogen_res_no = get_wo_intogen_res_no(experiment_no, experiment_df)    

        if np.isnan(intogen_res_no):
            res = res.assign(woIntogenPred=np.NaN)    
            res = res.assign(woIntogenPred_Rank=np.NaN)
        else: 
            intogen_res = pd.read_csv(res_dir + '/result_' + str(intogen_res_no) + '.tsv', sep='\t', index_col=0)
            intogen_res = aggregate_res_by_gene(intogen_res)
            intogen_res = intogen_res.rename(columns={'Prediction':'woIntogenPred'})
            res = pd.merge(res, intogen_res[['GeneID', 'woIntogenPred']], on='GeneID')
            res = res.assign(woIntogenPred_Rank = res.woIntogenPred.rank(ascending=False, method='min'))
    else:
        res = res.assign(woIntogenPred=np.NaN)
        res = res.assign(woIntogenPred_Rank=np.NaN)
        
    return res


def add_wo_outrider_res(res, outlier_input_feature, experiment_no, experiment_df, res_dir):
    if outlier_input_feature != '':
        outrider_res_no = get_wo_outrider_res_no(experiment_no, experiment_df)    

        if np.isnan(outrider_res_no):
            res = res.assign(woOutlierPred=np.NaN)    
            res = res.assign(woOutlierPred_Rank=np.NaN)
        else: 
            outrider_res = pd.read_csv(res_dir + '/result_' + str(outrider_res_no) + '.tsv', sep='\t', index_col=0)
            outrider_res = aggregate_res_by_gene(outrider_res)
            outrider_res = outrider_res.rename(columns={'Prediction':'woOutlierPred'})
            res = pd.merge(res, outrider_res[['GeneID', 'woOutlierPred']], on='GeneID')
            res = res.assign(woOutlierPred_Rank = res.woOutlierPred.rank(ascending=False, method='min'))
    else:
        res = res.assign(woOutlierPred=np.NaN)
        res = res.assign(woOutlierPred_Rank=np.NaN)
        
    return res


def add_wo_coess_res(res, coess_input_feature, experiment_no, experiment_df, res_dir):
    if coess_input_feature != '':
        coess_res_no = get_wo_coess_res_no(experiment_no, experiment_df)    

        if np.isnan(coess_res_no):
            res = res.assign(woCoessPred=np.NaN)    
            res = res.assign(woCoessPred_Rank=np.NaN)
        else: 
            coess_res = pd.read_csv(res_dir + '/result_' + str(coess_res_no) + '.tsv', sep='\t', index_col=0)
            coess_res = aggregate_res_by_gene(coess_res)
            coess_res = coess_res.rename(columns={'Prediction':'woCoessPred'})
            res = pd.merge(res, coess_res[['GeneID', 'woCoessPred']], on='GeneID')
            res = res.assign(woCoessPred_Rank = res.woCoessPred.rank(ascending=False, method='min'))
    else:
        res = res.assign(woCoessPred=np.NaN)
        res = res.assign(woCoessPred_Rank=np.NaN)
        
    return res


def get_features_post(features_full, sample_group, intogen_input_feature, outlier_input_feature):
    
    features_post = get_feature_model(features_full, sample_group, intogen_input_feature, outlier_input_feature, "")
    features_post = pd.concat([features_full[['gene_id']], features_post], axis=1)
    features_post = features_post.rename(columns={'gene_id': 'GeneID'})
    
    return features_post


def add_prediction_rank_post(res, intogen_input_feature, outlier_input_feature, coess_input_feature, features_post):
    if ((outlier_input_feature != '') | (intogen_input_feature != '')) * (coess_input_feature != ''):
        coess_select = features_post.apply(lambda x:sum(x[1:]), axis=1) == 0
        coess_select = pd.concat([features_post.GeneID, coess_select], axis=1).rename(columns={0:'only_coess'})
        # coess_select

        res = pd.merge(res, coess_select, on='GeneID', how='left')
        res = res.assign(Prediction_post=res.Prediction)
        res.loc[res.only_coess==True, 'Prediction_post'] = 0
        res = res.sort_values(by='Prediction_post', ascending=False)
        res = res.assign(Rank_post = res.Prediction_post.rank(ascending=False, method='min'))
    else:
        res = res.assign(Prediction_post=res.Prediction)
        res = res.assign(Rank_post = res.Prediction_post.rank(ascending=False, method='min'))
        res = res.assign(only_coess = np.NaN)
        
    return res


def calculate_rank_proportion(res_info):
    # calculate proportion
    rp_df = pd.DataFrame(columns=['Rank', 'Criteria', 'Proportion', 'Processing'])

    # create proportion df for rank
    rp_sub_pre_df = pd.DataFrame(columns=['isLabel', 'isCGC', 'isLeukemia', 'isOCG', 'isTSG', 'isIntOGen'], index=res_info.Rank.unique())

    # calculate proportion for rank
    for i in res_info.Rank.unique():
        rp_sub_pre_df.loc[i, 'isLabel'] = float(res_info[res_info.Rank <= i][['Label']].sum()/i)
        rp_sub_pre_df.loc[i, 'isCGC'] = float(res_info[res_info.Rank <= i][['isCGC']].sum()/i)
        rp_sub_pre_df.loc[i, 'isLeukemia'] = float(res_info[res_info.Rank <= i][['isLeukemia']].sum()/i)
        rp_sub_pre_df.loc[i, 'isOCG'] = float(res_info[res_info.Rank <= i][['isOCG']].sum()/i)
        rp_sub_pre_df.loc[i, 'isTSG'] = float(res_info[res_info.Rank <= i][['isTSG']].sum()/i)
        rp_sub_pre_df.loc[i, 'isIntOGen'] = float(res_info[res_info.Rank <= i][['isIntOGen']].sum()/i)
        
    # melt proportion for plotting for rank
    rp_sub_pre_df = rp_sub_pre_df.reset_index().rename(columns={'index':'Rank'})
    rp_sub_pre_df = pd.melt(rp_sub_pre_df, id_vars=['Rank'], var_name='Criteria', value_name='Proportion')
    rp_sub_pre_df['Criteria'] = pd.Categorical(rp_sub_pre_df['Criteria'], 
                                           categories=['isLabel', 'isCGC', 'isLeukemia', 'isOCG', 'isTSG', 'isIntOGen'], 
                                           ordered=True)
    rp_sub_pre_df['Processing'] = 'pre'

    rp_df = pd.concat([rp_df, rp_sub_pre_df])

    if all(res_info.Rank == res_info.Rank_post):
        rp_sub_post_df = rp_sub_pre_df
        rp_sub_post_df['Processing'] = 'post'
        rp_df = pd.concat([rp_df, rp_sub_post_df])
    else:    
        # create proportion df for rank_post
        rp_sub_post_df = pd.DataFrame(columns=['isLabel', 'isCGC', 'isLeukemia', 'isOCG', 'isTSG', 'isIntOGen'], index=res_info.Rank_post.unique())

        # calculate proportion for rank_post
        for i in res_info.Rank_post.unique():
            rp_sub_post_df.loc[i, 'isLabel'] = float(res_info[res_info.Rank_post <= i][['Label']].sum()/i)
            rp_sub_post_df.loc[i, 'isCGC'] = float(res_info[res_info.Rank_post <= i][['isCGC']].sum()/i)
            rp_sub_post_df.loc[i, 'isLeukemia'] = float(res_info[res_info.Rank_post <= i][['isLeukemia']].sum()/i)
            rp_sub_post_df.loc[i, 'isOCG'] = float(res_info[res_info.Rank_post <= i][['isOCG']].sum()/i)
            rp_sub_post_df.loc[i, 'isTSG'] = float(res_info[res_info.Rank_post <= i][['isTSG']].sum()/i)
            rp_sub_post_df.loc[i, 'isIntOGen'] = float(res_info[res_info.Rank_post <= i][['isIntOGen']].sum()/i)

        # melt proportion for plotting for rank_post
        rp_sub_post_df = rp_sub_post_df.reset_index().rename(columns={'index':'Rank'})
        rp_sub_post_df = pd.melt(rp_sub_post_df, id_vars=['Rank'], var_name='Criteria', value_name='Proportion')
        rp_sub_post_df['Criteria'] = pd.Categorical(rp_sub_post_df['Criteria'], 
                                               categories=['isLabel', 'isCGC', 'isLeukemia', 'isOCG', 'isTSG', 'isIntOGen'],
                                               ordered=True)
        rp_sub_post_df['Processing'] = 'post'
        rp_df = pd.concat([rp_df, rp_sub_post_df])
    
    return rp_df

def add_cancer_driver_info_feat(feat, cgc_cancer_gene, cgc_leukemia_gene, cgc_AML_gene, intogen_cancer_gene):
    feat = feat.assign(isCGC = np.isin(feat['gene_id'], cgc_cancer_gene['ENSGid']))
    feat = feat.assign(isLeukemia = np.isin(feat['gene_id'], cgc_leukemia_gene['ENSGid']))
    feat = feat.assign(isAML = np.isin(feat['gene_id'], cgc_AML_gene['ENSGid']))
    feat = feat.assign(isOCG = np.isin(feat['gene_id'], cgc_cancer_gene.ENSGid[cgc_cancer_gene.RoleinCancer.str.contains('oncogene').fillna(False)]))
    feat = feat.assign(isTSG = np.isin(feat['gene_id'], cgc_cancer_gene.ENSGid[cgc_cancer_gene.RoleinCancer.str.contains('TSG').fillna(False)]))
    feat = feat.assign(isIntOGen = np.isin(feat['geneSymbol'], intogen_cancer_gene['GeneSymbol']))
    return feat

def calculate_ratio_coess(cluster_to_curate, features_coess, feature):
    features_coess_curate = features_coess[['geneSymbol', cluster_to_curate, feature]]
    cluster = features_coess_curate.loc[features_coess_curate[cluster_to_curate] == 1]
    ratio = sum(cluster[feature])/len(cluster[feature])
    return ratio

def calculate_size_coess(cluster_to_curate, features_coess):
    features_coess_curate = features_coess[['geneSymbol', cluster_to_curate]]
    cluster = features_coess_curate.loc[features_coess_curate[cluster_to_curate] == 1]
    return len(cluster)

def calculate_coess_curate(coess_curate, features_coess):
    coess_curate['size'] = coess_curate.apply(lambda x: calculate_size_coess(x['feature'], features_coess), axis=1)
    coess_curate['label_ratio'] = coess_curate.apply(lambda x: calculate_ratio_coess(x['feature'], features_coess, 'isLabel'), axis=1)
    coess_curate['cgc_ratio'] = coess_curate.apply(lambda x: calculate_ratio_coess(x['feature'], features_coess, 'isCGC'), axis=1)
    coess_curate['leu_ratio'] = coess_curate.apply(lambda x: calculate_ratio_coess(x['feature'], features_coess, 'isLeukemia'), axis=1)
    coess_curate['ocg_ratio'] = coess_curate.apply(lambda x: calculate_ratio_coess(x['feature'], features_coess, 'isOCG'), axis=1)
    coess_curate['tsg_ratio'] = coess_curate.apply(lambda x: calculate_ratio_coess(x['feature'], features_coess, 'isTSG'), axis=1)
    return coess_curate

