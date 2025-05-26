# iLR </br> 
## 1. Introduction  
  iLR is a method selecting small sets of informative genes that distinguish subtle differences between cell states (i.e., disease, treatment) by iteratively selecting genes using a logistic regression framework coupled with penalized Pareto front optimization.
 

  
## 2. Installation
Depends: 

       Python (>= 3.9.2)

Requirements: 

      scanpy >= 1.8.2
      numpy >= 1.21.6
      pandas >= 1.5.3
      sklearn >= 1.5.3
      
  
## 3. Quick start

Run `iLR.py`. The parameters can be changed as below.

### 3.1 Prepare data

Anndata object can be accepted directly as input. `obs_names` corresponds to cells, and `var_names` corresponds to genes. Lognorm transformed expression and preselection of genes by the Wilcoxon test are recommended. 

Example of preparing data (from paper: https://www.science.org/doi/10.1126/science.aav8130, L4 cell type only):

     info = pd.DataFrame(columns = ['train_ASD_number', 'train_control_number', 'num_genes_started','test_ASD_number', 'test_control_number'])
     adata = L4

     sc.tl.rank_genes_groups(adata, groupby="diagnosis", method = 'wilcoxon')
     df = sc.get.rank_genes_groups_df(adata, group="ASD")
     selected = df[df['pvals_adj'] < 0.05]['names'].tolist()
     info.at[0,'num_genes_started'] = len(selected)

     X_all = adata.X
     y_all = np.asarray(adata.obs["diagnosis"])
     rs = ShuffleSplit(n_splits=1, test_size=0.3, random_state=0)
     rs.get_n_splits(X_all, y_all)
     train_index, test_index = next(rs.split(X_all, y_all)) 
     train_adata = adata[train_index,:]
     test_adata = adata[test_index,:]

     info.at[0,'train_ASD_number'] = np.sum(train_adata.obs['diagnosis'] == 'ASD')
     info.at[0,'train_control_number'] = np.sum(train_adata.obs['diagnosis'] == 'Control')
     info.at[0,'test_ASD_number'] = np.sum(test_adata.obs['diagnosis'] == 'ASD')
     info.at[0,'test_control_number'] = np.sum(test_adata.obs['diagnosis'] == 'Control')

     info
     
    |   train_ASD_number |   train_control_number |   num_genes_started |   test_ASD_number |   test_control_number |
    |--------------------|------------------------|---------------------|-------------------|-----------------------|
    |               2531 |                   2250 |                2941 |              1057 |                   993 |

 
Preprocessed demo data available at ...

### 3.2 Run scIAE
`iLR` returns a dictionary with penalty as key and corresponding dictionary of genes as values. Its inputs are listed below.

  - `train_data`: Anndata object with pre-filtered genes with lognorm transformed .X
  - `test_data`: Anndata object with pre-filtered genes with lognorm transformed .X
  - `observation`: obs name of interest (i.e., treatments)
  - `num_repeats`: number of repeats of iLR (default: 1)
  - `ia`: L2 regularization parameter (default: 0.1)
  - `per_remove`: proportion of genes removed at each iteration of iLR (default: 0.2)
  - `e`: penalty of Pareto Front (list, default: [0])
  - `min_num`: minimal number of genes desired in the final gene set (default: 10)
  - `rs`: random seed of logistic regression
  - `plot`: if plot the gene number vs classification AUC scatter plot (boolean, default: True)

 
 Run `iLR()`, then it will output evaluation table containing AUCs and gene sets with different pareto front penalty.
      
      counts, ev  = iLR(adata_train_filtered, adata_test_filtered, 'diagnosis',  ia = 0.1, e = [0, 1, 2], min_num = 10, plot = False)
      ev
     |   penalty |   train_mean_cv_acc |   train_mean_cv_auc |   number_selected_genes |   train_acc |   train_auc |   test_acc |   test_auc |
     |-----------|---------------------|---------------------|-------------------------|---    ----------|-------------|------------|------------|
     |         0 |            0.909224 |            0.972659 |                     314 |    0.942062 |    0.986966 |   0.830732 |   0.911936 |
     |         1 |            0.828067 |            0.908705 |                      81 |    0.8379   |    0.917084 |   0.800976 |   0.881176 |
     |         2 |            0.752351 |            0.831734 |                      25 |    0.756327 |    0.835111 |   0.745366 |   0.825654 |


Evaluation table contains:
   
   - `penalty`: penalty of Pareto front
   - `train_mean_cv_acc/auc`: 5-CV accuracy or AUC at the selected gene set on training dataset
   - `number_selected_genes`: number of genes selected
   - `train_acc/auc`: accuracy or AUC on training dataset
   - `test_acc/auc`: accuracy or AUC on testing dataset

    counts[2]
    
    Counter({'SPDYE2': 1,
         'GTF2H2': 1,
         'FAM73B': 1,
         'PELI2': 1,
         'PEBP1': 1,
         'RP11-711K1.8': 1,
         'TRABD2A': 1,
         'RP11-611E13.2': 1,
         'CLTCL1': 1,
         'NDUFAB1': 1,
         'BCYRN1': 1,
         'FAM153A': 1,
         'FBLN7': 1,
         'TRIML2': 1,
         'OR2L13': 1,
         'HSPA1A': 1,
         'LINC01482': 1,
         'HNRNPA3P6': 1,
         'TNC': 1,
         'HSPB1': 1,
         'AC105402.4': 1,
         'MX2': 1,
         'RP11-159J3.1': 1,
         'NPAS4': 1,
         'ZNF208': 1})
         
      
If `num_repeats` > 1, the values in each sub-counter will show how many times of a gene appear out of number of repeats. Evaluation table will have results for each repeat. 


