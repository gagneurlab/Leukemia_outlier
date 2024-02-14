import os
import sys
import time
import pandas as pd
import numpy as np
import xgboost as xgb

from functools import reduce
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV
from sklearn.utils.fixes import loguniform


def get_intogen_features(features_full, intogen_input_features, sample_group):
    if not bool(intogen_input_features):
        features_intogen = pd.DataFrame(index=features_full.index)
        return(features_intogen)
    pattern = "|".join(['^' + i + '-' for i in intogen_input_features.split(",")])
    features_intogen = features_full.filter(regex=pattern, axis=1).filter(like=sample_group, axis=1)
    return(features_intogen)


def get_outrider_features(features_full, outlier_input_features, sample_group):
    if not bool(outlier_input_features):
        features_outrider = pd.DataFrame(index=features_full.index)
        return(features_outrider)
    pattern = "|".join(['^' + i + '-' for i in outlier_input_features.split(",")])
    features_outrider = features_full.filter(regex=pattern, axis=1).filter(like=sample_group, axis=1)
    return(features_outrider)


def get_coess_features(features_full, coess_input_features):
    if not bool(coess_input_features):
        features_coess = pd.DataFrame(index=features_full.index)
        return(features_coess)
    pattern = "|".join(['^' + i + '-' for i in coess_input_features.split(",")])
    features_coess = features_full.filter(regex=pattern, axis=1)
    return(features_coess)


def get_feature_model(features_full, sample_group, intogen_input_features, outlier_input_features, coess_input_features):
    features_intogen = get_intogen_features(features_full, intogen_input_features, sample_group)
    features_outrider = get_outrider_features(features_full, outlier_input_features, sample_group)
    features_coess = get_coess_features(features_full, coess_input_features)
    features_model = pd.concat([features_intogen, features_outrider, features_coess], axis=1)
    return(features_model)


def print_res_path(res_full_path):
    print("")
    print('Result path:')
    print(res_full_path)

def initialize_res_df():
    res_df = pd.DataFrame(columns=['GeneID', 'GeneSymbol', 'Label', 'Prediction', 'random_repeat'])
    return res_df


def print_coeff_path(coeff_full_path):
    print("")
    print('Coefficiant path:')
    print(coeff_full_path)


def intitialize_coeff_df(model_method):
    if model_method == "lr":
        coeff_df = pd.DataFrame(columns=['Feature', 'Coeff', 'random_repeat'])
    if model_method == "rf":
        coeff_df = pd.DataFrame(columns=['Feature', 'Imp_mean_decrease', 'Imp_perm', 'random_repeat'])
    if model_method == "xgb":
        coeff_df = pd.DataFrame(columns=['Feature', 'feature_importances'])
    if model_method == "xgb_op":
        coeff_df = pd.DataFrame(columns=['Feature', 'feature_importances'])
    if model_method == "nn":
        coeff_df = pd.DataFrame(columns=['Feature', 'feature_importances'])    
    return(coeff_df)


def read_in_label_gene(label_dir_path, label_gene_list):
    driver_gene_list_path = os.path.join(label_dir_path, label_gene_list + "_list.tsv")
    print("")
    print('Label genes (cancer driver) path:')
    print(driver_gene_list_path)    
    driver_gene_list = pd.read_csv(driver_gene_list_path, sep='\t')
    driver_genes = driver_gene_list['ENSGid'].values
    return(driver_genes)


def create_label(features_full, driver_genes):
    features_genes = features_full['gene_id'].values
    labels = np.isin(features_genes, driver_genes)
    print("")
    print("Number of label genes/all genes in gencode v33b:")
    print(str(sum(np.isin(driver_genes, features_genes))) + 
          "/" + 
          str(len(features_genes)))
    return(labels)


def print_feature(feature):
    print("")
    print('Features:')
    print(feature.head())
    print("")
    print('Features dimension:')
    print(feature.shape)


def start_training():
    print("")
    print("")
    print("----- Training started -----")
    print("")
    print("")
    start_time = time.time()
    return(start_time)


def end_training(start_time):
    end_time = time.time()
    print("")
    print("--- Running time in total %s mins ---" % ( (end_time - start_time)/60 ))


def start_round(random_round):
    print("")
    print("")
    print('----- Random round:' + str(random_round) + ' -----')


def intialize_prediction_array(features):
    predictions = np.empty(len(features), dtype=float) 
    return(predictions)


def get_train_test_labels(labels, train_index_num, test_index_num):
    index_num = range(max(max(train_index_num), max(test_index_num)) + 1)
    train_index_logical = np.isin(index_num, train_index_num)
    test_index_logical = np.isin(index_num, test_index_num)
    train_labels, test_labels = labels[train_index_logical], labels[test_index_logical]
    
    print("")
    print('Train Features Index:', train_index_logical)  
    print('Num TRUE in Train Labels:', sum(train_labels))
    print('Test Features Index:', test_index_logical)
    print('Num TRUE in Test Labels:', sum(test_labels))
    return(train_labels, test_labels)


def get_train_test_features(features, train_index_num, test_index_num):
    index_num = range(max(max(train_index_num), max(test_index_num)) + 1)
    train_index_logical = np.isin(index_num, train_index_num)
    test_index_logical = np.isin(index_num, test_index_num)
    train_features, test_features = features.iloc[train_index_logical], features.iloc[test_index_logical]
    return(train_features, test_features)


def run_logistic_regression(features_model, labels,
                            train_index_num, test_index_num,
                            solver, penalty, max_iter, random_round):

    # get train/test
    train_labels, test_labels = get_train_test_labels(labels, train_index_num, test_index_num)
    train_features, test_features = get_train_test_features(features_model, train_index_num, test_index_num)
        
    model = LogisticRegression(solver=solver, penalty=penalty, max_iter=max_iter, random_state = random_round)
    clf = model.fit(train_features, train_labels)
    predictions_test = clf.predict_proba(test_features)[:,1]

    ### show attributes
    coeff_df_sub = pd.DataFrame({"Feature":features_model.columns, 
                                 "Coeff":clf.densify().coef_.flatten()}) 

    coeff_df_sub = pd.concat([
        coeff_df_sub, 
        pd.DataFrame({'Feature': 'Intercept','Coeff': clf.intercept_[0]}, index=[0])
    ], axis=0)
                             
    print('Coefficient:')
    print(coeff_df_sub)

    return(predictions_test, coeff_df_sub)


def run_random_forest(features_model, labels,
                      train_index_num, test_index_num,
                      n_estimators, min_samples_split, max_depth,
                      weight_true, weight_false, random_round):
    
    # get train/test
    train_labels, test_labels = get_train_test_labels(labels, train_index_num, test_index_num)
    train_features, test_features = get_train_test_features(features_model, train_index_num, test_index_num)
    
    rf = RandomForestClassifier(n_estimators = n_estimators, 
                                random_state = random_round,
                                min_samples_split = min_samples_split, 
                                max_depth = max_depth,
                                class_weight = {True: weight_true, False: weight_false})
    rf.fit(train_features, train_labels)
    predictions_test = rf.predict_proba(test_features)[:,1]

    coeff_df_sub = pd.DataFrame({"Feature": train_features.columns, 
                                 "Imp_mean_decrease": rf.feature_importances_,
#                                              "Imp_perm":result.importances_mean
                                   })
    print('Coefficient:')
    print(coeff_df_sub)
     
    return(predictions_test, coeff_df_sub)


def run_xgboost(features_model, labels,
                train_index_num, test_index_num,
                tree_method, random_round):
    
    # get train/test
    train_labels, test_labels = get_train_test_labels(labels, train_index_num, test_index_num)
    train_features, test_features = get_train_test_features(features_model, train_index_num, test_index_num)
    
    xgb_mod = xgb.XGBClassifier(tree_method=tree_method, random_state = random_round)
    xgb_mod.fit(train_features, train_labels)
    predictions_test = xgb_mod.predict_proba(test_features)[:,1]

    coeff_df_sub = pd.DataFrame({"Feature": train_features.columns, 
                                 "feature_importances": xgb_mod.feature_importances_
                                   })
    print('Coefficient:')
    print(coeff_df_sub)
     
    return(predictions_test, coeff_df_sub)


def run_xgboost_optimized(features_model, labels,
                          train_index_num, test_index_num,
                          tree_method, random_round):

    # get train/test
    train_labels, test_labels = get_train_test_labels(labels, train_index_num, test_index_num)
    train_features, test_features = get_train_test_features(features_model, train_index_num, test_index_num)
    
    param_dist = {
    "reg_lambda": loguniform(1e-2, 1e5), #define the search space for alpha and lambda
    "reg_alpha": loguniform(1e-2, 1e5)
    }

    cla = xgb.XGBClassifier(tree_method=tree_method, random_state = random_round)
    random_search = RandomizedSearchCV(
            cla, param_distributions=param_dist, n_iter=50, refit = True
    )
    random_search.fit(train_features, train_labels)
    predictions_test = random_search.predict_proba(test_features)[:,1]

    coeff_df_sub = pd.DataFrame({"Feature": train_features.columns, 
                                 "feature_importances": random_search.best_estimator_.feature_importances_
                                   })
    print('Coefficient:')
    print(coeff_df_sub)
     
    return(predictions_test, coeff_df_sub)


def run_nn(features_model, labels,
                          train_index_num, test_index_num, random_round, 
                          epochs, lr, input_dim, hidden_dim,
                          ):
    
    from torch.utils.data import Dataset, DataLoader
    from nn_model import nn_model, driver_dataset
    import torch
    
    # get train/test
    train_labels, test_labels = get_train_test_labels(labels, train_index_num, test_index_num)
    train_features, test_features = get_train_test_features(features_model, train_index_num, test_index_num)

    # Create data loaders
    train_dataset = driver_dataset(train_features, train_labels)
    test_dataset = driver_dataset(test_features, test_labels)
    
    train_loader = DataLoader(dataset=train_dataset, batch_size = len(train_features), shuffle=True)
    test_loader = DataLoader(dataset=test_dataset, batch_size = len(test_features), shuffle=False)
    
    model = nn_model(input_dim, hidden_dim, lr)
    model.set_writer(random_round, hidden_dim, lr, test_index_num)
    model.reset_weights(model)

    model.fit(train_loader, test_loader, epochs=epochs)
    predictions_test = model(torch.tensor(test_features.values, dtype = torch.float32))
    
    # no coefficients in NN, we only return the predictions
    return predictions_test.squeeze().detach().numpy()


def generate_result_table(features_full, labels, predictions, random_round):
    res_df = features_full[['gene_id', 'geneSymbol']]
    res_df.columns = ['GeneID', 'GeneSymbol']

    res_df = res_df.assign(Label = labels)
    res_df = res_df.assign(Prediction = predictions)

    res_df = res_df.sort_values(by=['Prediction'], ascending=False)
    res_df['random_repeat'] = random_round 

    print("")
    print(res_df.head())
    return res_df


def get_param_training(experiment_no, experiment_df):
    
    experiment_df_sub = experiment_df.loc[experiment_df.experiment_no == int(experiment_no)]
    
    label_gene_list =  experiment_df_sub['label_gene_list'].values[0]
    sample_group =  experiment_df_sub['sample_group'].values[0]
    random_seeds =  experiment_df_sub['random_seeds'].values[0]
    model_method =  experiment_df_sub['model_method'].values[0]
    
    return((label_gene_list, sample_group, random_seeds, model_method))


def get_param_feature(experiment_no, experiment_df):
    
    experiment_df_sub = experiment_df.loc[experiment_df.experiment_no == int(experiment_no)]
    
    intogen_input_feature =  experiment_df_sub['intogen_input_feature'].values[0]
    outlier_input_feature =  experiment_df_sub['outlier_input_feature'].values[0]
    coess_input_feature =  experiment_df_sub['coess_input_feature'].values[0]

    return((intogen_input_feature, outlier_input_feature, coess_input_feature))


def get_param_lr(experiment_no, experiment_df):
    
    experiment_df_sub = experiment_df.loc[experiment_df.experiment_no == int(experiment_no)]
    
    solver =  experiment_df_sub['solver'].values[0]
    penalty =  experiment_df_sub['penalty'].values[0] 
    max_iter =  experiment_df_sub['max_iter'].values[0]

    return((solver, penalty, max_iter))


def get_param_rf(experiment_no, experiment_df):
    
    experiment_df_sub = experiment_df.loc[experiment_df.experiment_no == int(experiment_no)]
    
    n_estimators =  experiment_df_sub['n_estimators'].values[0]
    min_samples_split =  experiment_df_sub['min_samples_split'].values[0]
    max_depth =  experiment_df_sub['max_depth'].values[0]

    weight_true =  experiment_df_sub['weight_true'].values[0]
    weight_false =  experiment_df_sub['weight_false'].values[0]

    return((n_estimators, min_samples_split, max_depth, weight_true, weight_false))


def get_param_xgb(experiment_no, experiment_df):
    
    experiment_df_sub = experiment_df.loc[experiment_df.experiment_no == int(experiment_no)]
    
    tree_method =  experiment_df_sub['tree_method'].values[0]

    return((tree_method))