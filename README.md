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
1. The Akaike Information Criterion (AIC) and confusion matrix were evaluation tools in logistic regression. Both models had the same confusion matrix and similar results, so changing the convergence tolerance did not greatly affect model prediction. Only 5 out of 204 were misclassified (2 false positives and 3 false negatives)

2. Random Forest visualization, OOB estimate of error rate, and confusion matrix were used as evaluation tools in Random Forest. model_rf_3 with 50 trees had the smallest error rate of 5.39%, with a generally decreasing error rate as the number of trees deccreased. In this model, 11 out of 204 samples were misclassified (5 false positives and 6 false negatives)

3. Prediction probability, confusion matrix, and area under the curve were evaluation metrics in XGBoost. All three models showed similar performance, with 87.8% accuracy. 


**Feature Importance**
1. An odds ratio was used to understand which features increase or decrease the odds of the outcome: PC2 (increase), PC3 (decrease)
2. Importance scores were used to determine feature importance in Random Forest: PC2, PC3, and PC6 
3. The importance attribute of the xgboost package was used to determine feature importance in XGBoost based on gain, cover, and frequency: PC2, PC1, and PC3 in order of most gain, cover, and frequency. There was no difference in accuracy based on are under the curve and confusion matrix. 5 samples were misclassified with 4 false positives and 1 false negative.

It can be concluded that PC2 and PC3 have the most predictive power based on the odds ratio and improtance scores.

**Visualization and Mapping Genes**
![Variance Importance Plot](Figures/Variance_Importance_Plot.png)
The variance importance plot from the random forest model shows that PC2 and PC3 contribute the most to the tumor classification.
![Barplot of Pathways](Figures/Barplot.png)
![Dotplot of Pathways](Figures/Dotplot.png)
Based on the Kegg results, most of the differentially expressed genes in the tumors are linked to pathways of neurodegeneration, cadherin signaling, and Alzheimer's disease. In fact, neurodegeneration and cancer stem from inverse cell patterns, where the former is linked to premature apoptosis and the latter is cell proliferation. This inverse correlation suggests that patients with Alzheimer's or Parkinson's disease have decreased susceptibility to lung cnacer because of the opposite cell cycle regulation patterns. Both disease types rely on tumor suppressor genes like p53. However, they still share similarities like increased risk with age (40-80). Moreover, both disease types are started by abnormal mitogenic triggers that intervene in the cell cycle. However, the barplot and dotplot show the least gene association with DNA replication and platelet activation, despite many cancer genes upregulating DNA replication to increase cell division. Moreover, the original study conducted through this dataset discussed the thromboietic programming involved with lung cancer, where the general trend is hyperactivation of platelets. Thrombopoietic programming is blueprint for this activation, but this dataset's pathway analysis illustrates low platelet activation compared to other pathways. This is suprising because platelet hyperactivation is common in lung cancer as an immunosuppressant and driver of TGF-Beta 1, which promots cell proliferation, especially in cases of tissue injury nad wound healing. Increased platelets generally encourages metastasis and angiogenesis, where tumors are nourished by blood vessels and invade neighboring oxygen supply. Therefore, further RNAseq datasets must be explored to map more differentially expressed genes and re-evaluate potentially the same genes to make a more confident biological conclusion about the role of DNA replication and platelet activation specifically in lung cancer.


**Discussion**
More hyperparameters shoudl be explored to see how performance fluctuates. Overall, random forest has the highest accuracy, with model_rf_3 of 50 trees having an error rate of 5.39%. However, the modeling cycle stages must be revisited to improve preprocessing with other encoding methods, organize the pipeline more, and try to minimize false negatives. Trusting a model with false negatives, which all three had, is dangerous in terms of diagnosing a patient. Therefore, while these models need to be refined, retrained, and re-evaluated with more data, they cannot be trusted blindly for a diagnosis. The unexpected results in terms of DNA replication and platelet activation genes being low in number also trigger the need for further investigation and analysis using core data science principles that align with the data and modeling cycle. However, the overall accuracy of all the models were fairly high with over 87%. Therefore, paired with human critical thinking, it can be used as a resource to facilitate diagnosis and treatment plans as long as it is accompanied by human scientific testing.
