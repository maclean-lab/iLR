import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
from sklearn.utils import shuffle
from sklearn.metrics import roc_auc_score
import math   



def iLR(adata, test_adata, observation,
             num_repeats = 1, ia = 0.1, per_remove = 0.2, e = [0], min_num = 10,
             rs = 1, plot = True):
    sc.set_figure_params(figsize=(10, 6))

    intToCounter = {}
    for k in e:
        intToCounter[k] = Counter()
    
    ev = pd.DataFrame(columns = ['penalty', 'train_mean_cv_acc', 'train_mean_cv_auc', 'number_selected_genes', 'train_acc', 
                                'train_auc', 'test_acc', 'test_auc'])
    
    
    for i in range(num_repeats):
        adata_filtered = adata
    
        scores = pd.DataFrame() 
        X_all = adata_filtered.X
        y_all = np.asarray(adata_filtered.obs[observation])
        X_all, y_all = shuffle(X_all, y_all)
        gene_list = []
        acc = []
        auc = []
        n_gene = []
        removed = []
        
        gene_selected =  np.asarray(adata_filtered.var_names)
        
        for j in range(100):

            log_reg = LogisticRegression(penalty='l2', C=ia, solver='liblinear', random_state = rs)
            scores['L2 logistic'] = cross_val_score(log_reg, X_all, y_all, cv=5)
            cv_mean_acc = np.mean(scores['L2 logistic'].to_list())
            acc.append(cv_mean_acc)
            scores['AUC'] = cross_val_score(log_reg, X_all, y_all, cv=5, scoring = 'roc_auc')
            cv_mean_auc = np.mean(scores['AUC'].to_list())
            auc.append(cv_mean_auc)
            
            l_ind = X_all.shape[1]
            n_gene.append(l_ind)

            
            log_reg_fit = log_reg.fit(X_all, y_all)
            importance = np.ravel(log_reg_fit.coef_)
            
            # use 5CV mean accuracy as cutoff
            n = int(np.ceil(l_ind * per_remove))
            indices = np.argpartition(abs(importance), n)
            ind = indices[n:]
            ind_rest = indices[:n]
            gene_selected = np.asarray(adata_filtered.var_names)[ind]
            gene_removed = np.asarray(adata_filtered.var_names)[ind_rest]

            removed.append(gene_removed.tolist())
            
            # check if reach min gene number
            if len(gene_selected) < min_num:
                break
        
            adata_filtered = adata_filtered[:, gene_selected]
            X_all = adata_filtered.X
            y_all = np.asarray(adata_filtered.obs[observation])
            X_all, y_all = shuffle(X_all, y_all)
            
        max_auc = np.max(auc)
        min_ngene = np.min(n_gene)
        st_ngene = (n_gene - np.mean(n_gene))/np.std(n_gene)
        std_auc = np.std(auc)
        if std_auc == 0:
            std_auc = 1
        st_auc = (auc - np.mean(auc))/std_auc
        dist = [math.dist([st_ngene[-1], np.max(st_auc)], [x, y]) for x, y in zip(st_ngene, st_auc)]
        min_dist = np.min(dist)
        
        for p in e:
            idx_list = [i for i, j in enumerate(dist) if min_dist <= j <= min_dist+p]
            idx = idx_list[-1]
                         
            gene_list = list(np.concatenate(removed[idx:]))
            gene_list.extend(gene_selected)
            
            X_train = adata[:,gene_list].X
            y_train = np.asarray(adata[:,gene_list].obs[observation])
            X_test = test_adata[:,gene_list].X
            y_test = np.asarray(test_adata[:,gene_list].obs[observation])

            log_reg = LogisticRegression(penalty='l2', C=0.1, solver='liblinear', random_state = rs).fit(X_train, y_train)
            train_acc = log_reg.score(X_train, y_train)
            train_auc = roc_auc_score(y_train, log_reg.decision_function(X_train))
    
            test_acc = log_reg.score(X_test, y_test)
            test_auc = roc_auc_score(y_test, log_reg.decision_function(X_test))


            df = {'penalty': p, 
             'train_mean_cv_acc': acc[idx], 
             'train_mean_cv_auc': auc[idx], 
             'number_selected_genes': len(gene_list), 
             'train_acc' : train_acc, 'train_auc' : train_auc, 
             'test_acc': test_acc, 'test_auc': test_auc}

            ev = ev._append(df, ignore_index = True)

            intToCounter[p].update(gene_list)
            
    
    if plot == True:
        # Plot the data
        plt.scatter(n_gene, auc, label = '5-cross validation mean auc')
        plt.gca().invert_xaxis()
        plt.xlabel("number of genes")
        #plt.ylabel("5-cross validation accuracy")
        #plt.legend()
        plt.show()
   
    

    return intToCounter, ev


