# Lung Cancer Classification in R
This R machine learning pipeline uses logistic regression, random forests, and XGBoost to classify normal and malignant tumors of a lung cancer dataset (GSE329380) from the Gene Expression Omnibus (GEO) of NCBI using RNA-Seq data.  The pipeline ends with a pathway enrichment analysis to map genes.

# Background

Lung cancer is ranked #1 globally as the leading cancer that causes death. Causes include smoking, radon, abestos, and other chemical carcinogens that cause malignant epigenetic changes that intervene functional cell checkpoints. While not all tumors are malignant, understanding differential expression patterns in various genes helps scientists to target relevant pathways when juxtaposed with normal transcriptomic data. Treatment usually includes readiation, surgery (e.h. lobectomy, pneumonectomy, resection), chemotherapy, stereotactic body radiotherapy, immunotherapy, and more.

**Dataset Description**

The machine learning pipeline uses the bulk RNA-seq data from an equal proportion of tumor and non-tumor lung tissue samples with 61852 genes and 204 patient samples. The original study used this data to understand the role of thrombopoietic programs in lung cancer, which drive blood platelet production.  


**GEO LINK**
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329380

# Motivation 

Understanding the genetic differences between normal and fatal tumors and accurately classifying tumors based on past data is important to avoid test errors and risk patient lives. Especially as contemporary society moves towards computational modeling, diagnosis, and innovative medicine production, training reliable models to facilitate treatment plans and tumor characteristics can allow healthcare professionals to spend more time on harnessing their technical expertise to support patients. In fact, targeted treatment through immunotherapy, AI-drived antibody discovery, and gene therapy are current examples of technology-driven innovative medicine.
This largely relies on existing expression data and mapping genes to pathways that explain other fatalities.

# Pipeline

**Data Preprocessing**

Label encoding: Converted gene expression counts from a dataframe to a matrix to perform numeric operations (e.g. isolating genes that are expressed in more than 3 samples).

Checked for null values.

Defined labels and expression data as a dataframe for differential expression analysis 

One-Hot encoding was employed to define malignant tumors as "1" and normal tissues as "0".

**Differential Expression Analysis**
Analyzed genes that were overexpressed or underexpressed in tumor samples to understand what genes are involved with tumor proliferation.
This was done using the DESeq2 package from Bioconductor.

**Normalization and Dimensionality Reduction**
Raw numerical values were adjusted to accurately compare gene expression, removing technical noise while preserving true biological signals usign the counts per million function.
Principal component analysis (PCA) was performed to compress the data into components, the top 10 of which could be used as the new features for model training.

**Model Training**
1. Logistic Regrsesion
This is one of the interpretable classification models that fits a sigmoid function on the data. Based on a probability score and threshold (0.5), the sample is classified as malignant or normal. 

The epsilon was treated as a hyperparameter representative of convergence tolerance where a model's improvement slows down (1e-8 and 1e-10).

2. Random Forest is a form of ensemble learning in machine learning that combines decision trees to perform both regression and classificatino tasks, choosing the average prediction in the former task and the majority prediction in the latter task. A forest of trees are created through bagging or bootstrap aggregating, sampling with replacement and choosing random features for each subset of data. The majority vote for a tumor class was considered in the prediciton.

The number of trees was treated as a hyperparameter and stopping condition (500, 100, 50).

3. XGBoost or eXtreme Gradient Boosting is a gradient boosted decision tree algorithm used for classification and regression. With parallel processing, the model builds decision trees sequentially, such that the new tree learns from the mistakes of the previous tree through the residuals. Using a gradient algorithm, the aggregation of trees annd computation of residuals allowed the model to improve and find tumor classifications and features that coverges to a minimum point. 

The number of rounds (nrounds) were treated as a hyperparameter (100, 500, 1000)

 

**Performance Evaluation**
1. The Akaike Information Criterion (AIC) and confusion matrix were evaluation tools in logistic regression.

2. Random Forest visualization, OOB estimate of error rae, and confusion matrix were used as evaluation tools in Random Forest.

3. Prediction probability, confusion matrix, and area under the curve were evaluation metrics in XGBoost.


**Feature Importance**
1. An odds ratio was used to understand which features increase or decrease the odds of the outcome: PC2 (increase), PC3 (decrease)
2. Importance scores were used to determine feature importance in Random Forest: PC2, PC3, and PC6 
3. The importance attribute of the xgboost package was used to determine feature importance in XGBoost based on gain, cover, and frequency: PC2, PC1, and PC3 in order of most gain, cover, and frequency.

It can be concluded that PC2 and PC3 have the most predictive power based on the odds ratio and improtance scores.

**Visualization**
```
![Alt text](https://github.com/aarthi073/Lung_Cancer_Classification.git/Figures/Variance_Importance_Plot.png)
```
** Mapping Genes **

