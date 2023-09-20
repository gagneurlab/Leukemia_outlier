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
# import functions
if os.getcwd().split('/')[-1] == 'vale':
    sys.path.append('Scripts/function/') # path to use when run the script from snakemake pipeline
else: 
    sys.path.append('../function/') # path to use when run the script locally from jupyterlab

from analysis_function import *

# %%
# read in snakemake object
file = open("/s/project/vale/driver_prediction_published_202309/processed_data/snakemake/postprocessing.p",'rb')
snakemake = pickle.load(file)

# %%
experiment_path = snakemake.input.experimentDesign
res_dir = os.path.dirname(snakemake.output.result_post)
project_dir = snakemake.params.projectPath
intogen_dir = snakemake.params.intogenDir

cgc_cancer_gene_path = snakemake.params.cgc_cancer_gene_processed
intogen_cancer_gene_path = snakemake.params.intogen_cancer_gene

# %% [markdown]
# # Define plots

# %%
label_gene_list = 'CGC_leukemia_gene'
sample_group = 'leukemia_14group'
model_method = 'rf'
intogen_input_feature = ['clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv', '']
outlier_input_feature = ['absplice','fr', 'or,ac', 'or,ac,absplice,fr', 'or,ac,fr']
coess_input_feature = ['']
plot_sub_dir = 'test'

save_plot = True

exp_tag_common_diaplay = False 
exp_tag_individual_display_list = False

# matplotlib color cycle
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# %%
# add more colors
color_cycle_customized = color_cycle + ["crimson", "indigo", "olive"]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_cycle_customized)
color_cycle = color_cycle_customized

# %%
# global plot setup for matplotlib and ggplot
plt.rcParams["figure.figsize"] = (8, 8)
font = {'size': 14}

plt.rc('font', **font)

# %% [markdown]
# # Preparation

# %%
plot_dir = project_dir + '/Output/plot/' + plot_sub_dir
if not os.path.isdir(plot_dir):
    os.makedirs(plot_dir)

experiment_design = pd.read_csv(experiment_path, sep='\t')
experiment_design = experiment_design.fillna("")

# %%
exp_viz = experiment_design[
    (np.isin(experiment_design.label_gene_list, label_gene_list)) *
    (np.isin(experiment_design.sample_group, sample_group)) *
    (np.isin(experiment_design.model_method, model_method)) *
    (np.isin(experiment_design.intogen_input_feature, intogen_input_feature)) *
    (np.isin(experiment_design.outlier_input_feature, outlier_input_feature)) *
    (np.isin(experiment_design.coess_input_feature, coess_input_feature))
]

exp_viz

# %%
exp_viz = simplify_exp(exp_viz)
exp_viz, exp_tag_common, exp_tag_common_diaplay, exp_tag_individual_display_list = add_exp_tag(exp_viz, exp_tag_common_diaplay, exp_tag_individual_display_list)
exp_viz = add_result_paths(exp_viz, project_dir, intogen_dir)

exp_viz

# %% [markdown]
# # Precision Recall Curve

# %%
# plot Recall and Precision
p_prc = plt.figure()

performance_df = pd.DataFrame(columns=['Training_setup', 'random_repeat', 'AP', 'auPRC']) 
prc_df = pd.DataFrame(columns=['Training_setup', 'random_repeat', 'precision','recall'])
prc_ci_df = pd.DataFrame(columns=['Training_setup', 'random_repeat', 'recall', 'precision_mean', 'precision_median', 'precision_std'])
    
for exp_viz_index, exp_viz_row in exp_viz.iterrows():
    # read in res
    res = pd.read_csv(exp_viz_row['res_pre_path'], sep='\t')
    exp_tag_individual = exp_viz_row['exp_tag_individual']
    exp_tag_individual_display = exp_viz_row['exp_tag_individual_display']
    
    performance_df_ts = pd.DataFrame(columns=['Training_setup', 'random_repeat', 'AP', 'auPRC']) 

    for rep in res.random_repeat.unique():
        
        # extract res of certain random_repeat
        res_sub = res[res['random_repeat'] == rep]
        labels = res_sub['Label']
        predictions = res_sub['Prediction']

        # calculate AP and auPRC
        average_precision = average_precision_score(labels, predictions)
        precision, recall, _ = precision_recall_curve(labels, predictions)
        area = auc(recall, precision)

        performance_df_rep = pd.DataFrame({'Training_setup': exp_tag_individual, 
                                           'random_repeat': rep,
                                           'AP': average_precision,
                                           'auPRC': area}, index=[0])
        performance_df_ts = pd.concat([performance_df_ts, performance_df_rep])
        performance_df = pd.concat([performance_df, performance_df_rep])
        
        prc_df_rep = pd.DataFrame({'Training_setup': exp_tag_individual, 
                                   'random_repeat': rep,
                                   'precision': precision,
                                   'recall': recall}) # unique for each random repeat
        prc_df = pd.concat([prc_df, prc_df_rep]) # unique for each training setup

        
    prc_df_ts = prc_df[prc_df.Training_setup == exp_tag_individual]
    average_precision_median = performance_df_ts.AP.median()
    
    for rc in prc_df_ts.recall.unique():
        
        # using all the precision per recall
        pr_per_rc = prc_df_ts[prc_df_ts.recall == rc].precision
        prc_ci_df_rc = pd.DataFrame({'Training_setup': exp_tag_individual, 
                                     'random_repeat': rep,
                                     'recall': rc, 
                                     'precision_mean': pr_per_rc.mean(),
                                     'precision_median': pr_per_rc.median(),
                                     'precision_std': pr_per_rc.quantile(),
                                     'precision_up': pr_per_rc.quantile(0.9),
                                     'precision_dn': pr_per_rc.quantile(0.1)}, index=[0])
        prc_ci_df = pd.concat([prc_ci_df, prc_ci_df_rc]) 

        
    # plot median 90%q 10%q
    prc_ci_df_ts = prc_ci_df[prc_ci_df.Training_setup == exp_tag_individual].sort_values(by='recall')
    prc_ci_df_ts = prc_ci_df_ts.astype({'recall':'float'})
    area = auc(prc_ci_df_ts['recall'], prc_ci_df_ts['precision_median'])
    
    if 'plot_label' in locals():
        plt.step(prc_ci_df_ts['recall'], 
                 prc_ci_df_ts['precision_median'], 
                 where='pre', 
                 label='AP median={0:0.4f}'.format(average_precision_median) + '  ' + plot_label[i])    
    else:
        plt.step(prc_ci_df_ts['recall'], 
                 prc_ci_df_ts['precision_median'], 
                 where='pre', 
                 label='AP median={0:0.4f}'.format(average_precision_median) + '  ' + exp_tag_individual_display)

    plt.fill_between(prc_ci_df_ts['recall'], prc_ci_df_ts['precision_dn'], prc_ci_df_ts['precision_up'], alpha=.1)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.0])
    plt.xlim([0.0, 1.0])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('Precision-Recall curve'+ exp_tag_common_diaplay, pad=15)
    plt.legend(loc="upper right")
    plt.grid()

plt.show()

# %%
if save_plot is True:
    p_prc_path = os.path.join(plot_dir, exp_tag_common + '-prc.png')
    p_prc.savefig(p_prc_path)

# %% [markdown]
# # Average Precision plot

# %%
# calculate median, 90% q, 10% q performance
ap_dfs =  [pd.pivot_table(performance_df, values='AP', index = ['Training_setup'], aggfunc=np.median).reset_index().rename(columns={'AP':'AP_median'}),
           pd.pivot_table(performance_df, values='AP', index = ['Training_setup'], aggfunc=lambda y: np.percentile(y, 90)).reset_index().rename(columns={'AP':'AP_9q'}),
           pd.pivot_table(performance_df, values='AP', index = ['Training_setup'], aggfunc=lambda y: np.percentile(y, 10)).reset_index().rename(columns={'AP':'AP_1q'})]
ap_df = reduce(lambda  left,right: pd.merge(left,right,on=['Training_setup']), ap_dfs)

# # add type informatioin
# ap_df['Type'] = ''
# ap_df.loc[['g7tools' in x for x in ap_df.Training_setup], 'Type'] = 'Genome'
# ap_df.loc[['or' in x for x in ap_df.Training_setup], 'Type'] = 'Transcriptome'
# ap_df.loc[['ac' in x for x in ap_df.Training_setup], 'Type'] = 'Transcriptome'
# ap_df.loc[['coess' in x for x in ap_df.Training_setup], 'Type'] = 'External'
# ap_df['Type'] = pd.Categorical(ap_df['Type'], categories=['Genome', 'Transcriptome', 'External'], ordered=True)

ap_df['Training_setup'] = pd.Categorical(ap_df['Training_setup'], 
                                         categories=exp_viz.exp_tag_individual.to_list(), 
                                         ordered=True)

# %%
p_ap = (ggplot(ap_df)
      + geom_line(aes(x='Training_setup', y='AP_median', group = 1), linetype='dashed')
      + geom_linerange(aes(x='Training_setup', ymax='AP_9q', ymin='AP_1q'))
      + geom_point(aes(x='Training_setup', y='AP_median', color = 'Training_setup', group = 1), size = 4, alpha = 0.6)
      
      + scale_color_manual(values=color_cycle)
      + scale_x_discrete(labels=exp_tag_individual_display_list)
      # + ylim(0, 0.2)
      
      + xlab('Training features')
      + ylab('Performance (AP)')
      + labs(color = "Training features")
      + ggtitle('Average Precision' + exp_tag_common_diaplay)
        
      + theme_bw()
      + theme(legend_position = "none")
      + theme(axis_text_x=element_text(size=10))
      + theme(axis_text_y=element_text(size=10))
      + theme(axis_text_x=element_text(rotation=45, hjust=1))
     )
                  
p_ap

# %%
if save_plot is True:
    p_ap_path = os.path.join(plot_dir, exp_tag_common + '-ap.png')
    p_ap.save(p_ap_path, height=4, width=6, units = 'in', dpi=300)

# %% [markdown]
# # Rank Proportion plot

# %%
rp_df = pd.DataFrame(columns=['Rank', 'Criteria', 'Proportion', 'Training_setup', 'Processing'])

for exp_viz_index, exp_viz_row in exp_viz.iterrows():
    # read in res
    rp_sub_df = pd.read_csv(exp_viz_row['rp_path'], sep='\t')
    rp_sub_df = rp_sub_df.assign(Training_setup = exp_viz_row['exp_tag_individual'])
    
    rp_df = pd.concat([rp_df, rp_sub_df])

rp_display = ''
pre_post_identical = np.array_equal(rp_df.loc[rp_df.Processing=='pre', 'Proportion'], rp_df.loc[rp_df.Processing=='post', 'Proportion'])
if pre_post_identical:
    rp_display = ' (pre/post identical)'


# %%
rp_plot_df = rp_df
# rp_plot_df = rp_plot_df.loc[rp_plot_df.Processing == 'post']
# rp_plot_df = rp_plot_df.loc[np.logical_or(rp_plot_df.Criteria == 'isCGC', rp_plot_df.Criteria == 'isIntOGen')]

rp_plot_df['Rank'] = rp_plot_df.Rank.astype(float)
rp_plot_df['Proportion'] = rp_plot_df.Proportion.astype(float)
rp_plot_df['Training_setup'] = pd.Categorical(rp_plot_df['Training_setup'], 
                                              categories=exp_viz.exp_tag_individual.to_list(), 
                                              ordered=True)
rp_plot_df['Criteria'] = pd.Categorical(rp_plot_df['Criteria'], 
                                        categories=['isLabel', 'isCGC', 'isLeukemia', 'isOCG', 'isTSG', 'isIntOGen'],
                                        ordered=True)

rp_plot_df

# %%
# to avoid error regarding numpy64
np.float = float  

# %% [markdown]
# ## all - log scale

# %%
p_rp_log10 = (ggplot(rp_plot_df, aes(x='Rank', y='Proportion', color='Training_setup'))
      + geom_line(size = 1, alpha=0.5)
#       + geom_step(direction='vh')
      + facet_grid('Processing ~ Criteria')
      
      + ylim(0,1)
      + scale_x_log10()
      + scale_color_manual(values=color_cycle)
              
      + ylab("Proportion of genes")
      + labs(color = "Training features")
      + ggtitle('Rank Proportion Plot' + exp_tag_common_diaplay + rp_display)
              
      + theme_bw()
      + theme(axis_text_x=element_text(size=10))
      + theme(axis_text_y=element_text(size=10))
      + theme(axis_text_x=element_text(rotation=45, hjust=1))
     )

if exp_tag_individual_display_list != False:
    p_rp_log10 = p_rp_log10 + scale_color_manual(labels=exp_tag_individual_display_list, values=color_cycle)
    
p_rp_log10

# %%
if save_plot is True:
    p_rp_log10_path = os.path.join(plot_dir, exp_tag_common + '-rp.png')
    p_rp_log10.save(p_rp_log10_path, height=6, width=12, units = 'in', dpi=300)

# %% [markdown]
# ## top200

# %%
n_to_plot = 200

p_rp_200 = (ggplot(rp_plot_df.loc[np.logical_and(rp_plot_df.Criteria == 'isCGC', rp_plot_df.Processing == 'post')], 
                   aes(x='Rank', y='Proportion', color='Training_setup'))
      + geom_line(size = 1, alpha=0.5)
#       + geom_step(direction='vh')
            
      + xlim(0, n_to_plot)
      + ylim(0,1)
      + scale_color_manual(values=color_cycle)
            
      + ylab("Proportion of CGC genes")
      + labs(color = "Training features")
      + ggtitle('Rank Proportion Plot' + exp_tag_common_diaplay + rp_display)
            
      + theme_bw()
      + theme(axis_text_x=element_text(size=10))
      + theme(axis_text_y=element_text(size=10))
      + theme(axis_text_x=element_text(rotation=45, hjust=1))
     )

if exp_tag_individual_display_list != False:
    p_rp_200 = p_rp_200 + scale_color_manual(labels=exp_tag_individual_display_list, values=color_cycle)
    
p_rp_200

# %%
n_to_plot = 200

p_rp_200 = (ggplot(rp_plot_df, aes(x='Rank', y='Proportion', color='Training_setup'))
      + geom_line(size = 1, alpha=0.5)
#       + geom_step(direction='vh')
      + facet_grid('Processing ~ Criteria')
      # + facet_wrap('Criteria', nrow = 1)
            
      + xlim(0, n_to_plot)
      + ylim(0,1)
      + scale_color_manual(values=color_cycle)
            
      + ylab("Proportion of genes")
      + labs(color = "Training features")
      + ggtitle('Rank Proportion Plot' + exp_tag_common_diaplay + rp_display)
            
      + theme_bw()
      + theme(axis_text_x=element_text(size=10))
      + theme(axis_text_y=element_text(size=10))
      + theme(axis_text_x=element_text(rotation=45, hjust=1))
     )

if exp_tag_individual_display_list != False:
    p_rp_200 = p_rp_200 + scale_color_manual(labels=exp_tag_individual_display_list, values=color_cycle)
    
p_rp_200

# %%
if save_plot is True:
    p_rp_200_path = os.path.join(plot_dir, exp_tag_common + '-rp'+str(n_to_plot)+'.png')
    p_rp_200.save(p_rp_200_path, height=6, width=12, units = 'in', dpi=300)

# %% [markdown]
# # area under Rank Propotion 200 plot

# %%
n_to_plot = 200
auRP_df = pd.DataFrame(columns = ['Rank', 'Criteria', 'Training_setup', 'Processing'])

iterables = [ rp_plot_df.Criteria.unique(), rp_plot_df.Training_setup.unique(), rp_plot_df.Processing.unique() ]

for t in itertools.product(*iterables):
    i = t[0]
    j = t[1]
    k = t[2]

    rp_plot_df_sub = rp_plot_df[
        (rp_plot_df.Criteria == i) *
        (rp_plot_df.Training_setup == j) *
        (rp_plot_df.Processing == k)]
    rp_plot_df_sub = rp_plot_df_sub[rp_plot_df_sub.Rank <= n_to_plot]
    rp_plot_df_sub = rp_plot_df_sub.sort_values('Rank')
    area = auc(rp_plot_df_sub.Rank, rp_plot_df_sub.Proportion)

    auRP_df_sub = pd.DataFrame({'Rank': [n_to_plot], 'Criteria': [i], 'Training_setup':[j], 'Processing': [k],
                                'auRP' + str(n_to_plot): [area]})
    auRP_df = pd.concat([auRP_df, auRP_df_sub]) 

auRP_df['averProp'+str(n_to_plot)] = auRP_df['auRP'+str(n_to_plot)]/n_to_plot
    
auRP_df['Training_setup'] = pd.Categorical(auRP_df['Training_setup'], 
                                           categories=exp_viz.exp_tag_individual.to_list(), 
                                           ordered=True)
auRP_df['Criteria'] = pd.Categorical(auRP_df['Criteria'], 
                                     categories=['isLabel', 'isCGC', 'isLeukemia', 'isOCG', 'isTSG', 'isIntOGen'],
                                     ordered=True)

# %% [markdown]
# ## pre/post 

# %%
p_auRP = (ggplot(auRP_df, aes(x='Training_setup', y='auRP'+str(n_to_plot), color = 'Training_setup'))
        + geom_point(size = 4, alpha = 0.6)
        + geom_line(aes(x='Training_setup', y='auRP'+str(n_to_plot), group = 1), linetype='dashed', color='black')
        + facet_grid('Processing ~ Criteria', scales='fixed')

        + scale_color_manual(values=color_cycle)

        + xlab('Training features')
        + ylab('Performance (auRP' + str(n_to_plot) + ')')
        + ggtitle('auRankProportion' + exp_tag_common_diaplay + rp_display)

        + theme_bw()
        + theme(legend_position = "none")
        + theme(axis_text_x=element_text(size=10))
        + theme(axis_text_y=element_text(size=10))
        + theme(axis_text_x=element_text(rotation=45, hjust=1))
     )

if exp_tag_individual_display_list != False:
    p_auRP = p_auRP + scale_x_discrete(labels=exp_tag_individual_display_list)

p_auRP

# %%
if save_plot is True:
    p_auRP_path = os.path.join(plot_dir, exp_tag_common + '-auRP'+str(n_to_plot)+'.png')
    p_auRP.save(p_auRP_path, height=6, width=15, units = 'in', dpi=300)

# %% [markdown]
# ## post

# %%
auRP_post_df = auRP_df.loc[auRP_df.Processing=='post']

# %%
p_auRP_post = (ggplot(auRP_post_df, aes(x='Training_setup', y='auRP'+str(n_to_plot), color = 'Training_setup'))
        + geom_point(size = 4, alpha = 0.6)
        + geom_line(aes(x='Training_setup', y='auRP'+str(n_to_plot), group = 1), linetype='dashed', color='black')
        + facet_wrap('Criteria', scales='fixed', nrow = 2)

        + scale_color_manual(values=color_cycle)

        + xlab('Training features')
        + ylab('Performance (auRP' + str(n_to_plot) + ')')
        + ggtitle('auRankProportion' + exp_tag_common_diaplay + rp_display)

        + theme_bw()
        + theme(legend_position = "none")
        + theme(axis_text_x=element_text(size=10))
        + theme(axis_text_y=element_text(size=10))
        + theme(axis_text_x=element_text(rotation=45, hjust=1))
     )

if exp_tag_individual_display_list != False:
    p_auRP_post = p_auRP_post + scale_x_discrete(labels=exp_tag_individual_display_list)

p_auRP_post

# %%
if save_plot is True:
    p_auRP_post_path = os.path.join(plot_dir, exp_tag_common + '-auRP'+str(n_to_plot)+'post.png')
    p_auRP_post.save(p_auRP_post_path, height=6, width=10, units = 'in', dpi=300)

# %% [markdown]
# ## post - Average Propotion

# %%
p_averProp_post = (ggplot(auRP_post_df, aes(x='Training_setup', y='averProp'+str(n_to_plot), color = 'Training_setup'))
        + geom_point(size = 4, alpha = 0.6)
        + geom_line(aes(x='Training_setup', y='averProp'+str(n_to_plot), group = 1), linetype='dashed', color='black')
        + facet_wrap('Criteria', scales='fixed', nrow = 2)

        + scale_color_manual(values=color_cycle)

        + xlab('Training features')
        + ylab('Performance (averProp' + str(n_to_plot) + ')')
        + ggtitle('Average Proportion within Rank ' + str(n_to_plot) + exp_tag_common_diaplay + rp_display)

        + theme_bw()
        + theme(legend_position = "none")
        + theme(axis_text_x=element_text(size=10))
        + theme(axis_text_y=element_text(size=10))
        + theme(axis_text_x=element_text(rotation=45, hjust=1))
     )

if exp_tag_individual_display_list != False:
    p_averProp_post = p_averProp_post + scale_x_discrete(labels=exp_tag_individual_display_list)

p_averProp_post

# %%
p_averProp_post = (ggplot(auRP_post_df.loc[auRP_post_df.Criteria=='isCGC'], 
                          aes(x='Training_setup', y='averProp'+str(n_to_plot), color = 'Training_setup'))
        + geom_point(size = 4, alpha = 0.6)
        + geom_line(aes(x='Training_setup', y='averProp'+str(n_to_plot), group = 1), linetype='dashed', color='black')

        + scale_color_manual(values=color_cycle)

        + xlab('Training features')
        + ylab('Performance (averProp' + str(n_to_plot) + ')')
        + ggtitle('Average Proportion within Rank ' + str(n_to_plot) + exp_tag_common_diaplay + rp_display)

        + theme_bw()
        + theme(legend_position = "none")
        + theme(axis_text_x=element_text(size=10))
        + theme(axis_text_y=element_text(size=10))
        + theme(axis_text_x=element_text(rotation=45, hjust=1))
     )

if exp_tag_individual_display_list != False:
    p_averProp_post = p_averProp_post + scale_x_discrete(labels=exp_tag_individual_display_list)

p_averProp_post

# %%
if save_plot is True:
    p_averProp_post_path = os.path.join(plot_dir, exp_tag_common + '-averProp'+str(n_to_plot)+'post.png')
    p_averProp_post.save(p_averProp_post_path, height=6, width=10, units = 'in', dpi=300)

