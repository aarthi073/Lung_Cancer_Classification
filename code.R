#Using RNASeq data to classify malignant tumors through logistic regression, a random forest, and XGBoost tree

#1
#install bioconductor (packages for gene experssion, genomics, RNA-seq, microarrays)
install.packages("BiocManager")
install.packages("randomForest")
install.packages("devtools")

BiocManager::install("GEOquery")
BiocManager::install("edgeR")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("ReactomePA")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq2")
devtools::install_github('araastat/reprtree')
reprtree::plot.getTree(model_rf, k = 1)
install.packages("xgboost")
install.packages("caret")
install.packages("pROC")


library(BiocManager)
library(GEOquery)
library(edgeR)
library(randomForest)
#look at tree
library(rpart)
library(reprtree)
library("xgboost")
library(caret)
library(pROC)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DESeq2)

#extract only the expression matrix of lung cancer dataset
gset <- getGEO("GSE329380", GSEMatrix = TRUE)
expr <- exprs(gset[[1]])
dim(expr)

# Load the lung cancer count matrix
counts <- read.csv("GSE329380/GSE329380_expr_count_lungcaner.csv.gz", header = TRUE)
gene_ids <- counts$gene_id
counts_mat <- as.matrix(counts[, -(1:2)]) #remove gene_id columns
storage.mode(counts_mat) <- "numeric"
rownames(counts_mat) <- gene_ids
# Inspect the dataset size: rows (genes) x columns (patient samples)
dim(counts_mat)

# Look at a slice of the raw integer reads
counts[1:5, 1:4]


#check for null values
clean <- any(is.na(counts))
dim(clean)
#2
#Differential Expression Analysis 

#Start at character 1 and stop at character 1; if malignant("C"), mark as 1; otherwise, mark as 0
sample_names <- colnames(counts_mat)

labels <- ifelse(substr(sample_names, 1, 1) == "C", "Cancer", "Normal")

table(labels)
dim(X)

colData <- data.frame(
  row.names = sample_names,
  condition = factor(labels, levels=c("Normal", "Cancer"))
)
tail(colData)
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = colData,
                              design = ~ condition)

#generate results based on DE dataset from previous lines 
dds <- DESeq(dds)


#confirm result names
resultsNames(dds)


res <- results(dds, contrast = c("condition", "Cancer", "Normal"))

#filter low-expression noise
#C for cancer, N for normal
#Expression of over 10 in at least 3 samples
genes <- rowSums(counts(dds) >= 10) >= 3
filtered_counts <- counts_mat[genes, ]
#3 
#NORMALIZATION -> adjust the raw numerical values so you can accurately compare gene expression,
#removing technical noise while preserving true biological signals.
#counts per million



normalized_matrix <- cpm(filtered_counts, log=TRUE, prior.count=1)
#PCA transposition
pca_input <- t(normalized_matrix)

#Run PCA
pca <- prcomp(pca_input, center=TRUE, scale.=TRUE)

summary(pca)

#4a Train Logistic Regression Model
#features: first 10 components
X <- pca$x[, 1:10]
#4a
#TRAINING LOGISTIC REGRESSION 
#numeric labels
labels <- ifelse(substr(rownames(pca_input), 1, 1) == "C", 1, 0)
table(labels)

#building dataset used in model
data <- data.frame(
  label = labels,
  X
)



 #logistic regression; Generalized Linear Model (linear regression to classification problems)
#binomial classification problem so#train family is binomial
#epsilon -> convergence tolerance (predictions stabilize); default is 1e-8
model_lgr_1 <- glm(label ~ ., data = data, family = binomial())


#Positive Coefficient increasing tumor probability, decreases tumor probability
summary(model_lgr_1)
print(model_lgr_1)
#predictions
probs_1 <- predict(model_lgr_1, type="response")
print(probs_1)

#convert probabilities
#if probability is greater than 0.5, it is a tumor. Otherwise, it is normal.
pred_1 <- ifelse(probs_1 > 0.5, 1, 0)

#5a evaluate performance; confusion matrix
table(Predicted = pred_1, Actual = labels)

#episilon = 1e-10; slightly lower AIC by 0.002 
model_lgr_2 <- glm(label ~ ., data = data, family = binomial(), epsilon = 1e-10)


summary(model_lgr_2)
print(model_lgr_2)
#predictions
probs_2 <- predict(model_lgr_2, type="response")
print(probs_2)

#convert probabilities
#if probability is greater than 0.5, it is a tumor. Otherwise, it is normal.
pred_2 <- ifelse(probs_2 > 0.5, 1, 0)

#5a evaluate performance; confusion matrix
table(Predicted = pred_2, Actual = labels)

#5b Random Forest Modeling
#can change number of trees as hyperparameter

#data_clean <- data[, !colnames(data) %in% "label"]


model_rf <- randomForest(factor(label) ~.,
                        data = data,
                        ntree = 500, 
                        mtry =3,
                        importance = TRUE)
summary(model_rf)


model_rf


new_model_rf <- randomForest(x=as.matrix(X),
                             y = factor(labels))

summary(new_model_rf)
    

new_model_rf

#100 trees
model_rf_2 <- randomForest(factor(label) ~.,
                         data = data,
                         ntree = 100, 
                         mtry =3,
                         importance = TRUE)
summary(model_rf_2)
model_rf_2

#50 trees
model_rf_3 <- randomForest(factor(label) ~.,
                           data = data,
                           ntree = 50, 
                           mtry =3,
                           importance = TRUE)
summary(model_rf_3)
model_rf_3



#3 Visualize tree
varImpPlot(model_rf)
#Random Forest Evaluation
#Tree without matrix has lower error

#5c XGBoost Modeling

#requires matrix input
X <- as.matrix(data[,-1])
y <- as.integer(as.character(data$label))

#split training and testing data
set.seed(42)
idx <- sample(nrow(data), 0.8 * nrow(data))
train_X <- X[idx, ]
train_y <- y[idx]
test_X <- X[-idx, ]
test_y <- y[-idx]
train <- xgb.DMatrix(data = train_X, label=train_y)
test <- xgb.DMatrix(data = test_X, label = test_y)

#hyperparameters: nrounds= 100, 500, 1000
model_xgb <- xgb.train(data=train,
                     nrounds = 100,
                     objective = "binary:logistic",
                     watchlist = list(train = train, test = test),
                     early_stopping_rounds = 50,
                     eval_metric = "auc",                    
                     verbose=0)

pred_prob <- predict(model_xgb, test) 
pred_class <- ifelse(pred_prob > 0.5, 1, 0)

#Evaluate XGBoost tree
confusionMatrix(factor(pred_class), factor(test_y))
auc(roc(test_y, pred_prob))


#6 feature importance logistic regression -> OR >1 means 
#Odds ratio
ivl <- exp(cbind(OR = coef(model_lgr_2), confint(model_lgr_2)))
#remove intercept
#ivl <- ivl[-1,]
ivl <- as.data.frame(exp(cbind(OR=coef(model_lgr_2), confint(model_lgr_2))))
print(ivl)
#OR>1 increases odds of outcome
sum(ivl$"2.5 %">1)
sum(ivl$"97.5 %">1)

#feature importance Random Forest
importance_scores <- importance(model_rf)
print(importance_scores)



#feature importance XGBoost
importance <- xgb.importance(feature_names = colnames(train_X), model=model_xgb)
print(importance)

#intercept, PC1, PC2, PC5, PC7, PC8  ???

#ensemble Gene IDs
#7 Map genes to pathways
#annotated databases
#Gene Ontology has a human gene annotation database and cluster profiler
#KEGG details metabolic and signaling pathways
#remove version number from gene_id by matching the literal period and anything after it
gene_ids_clean <- sub("\\..*", "", gene_ids)

#KEGG works with Entrez IDs (integer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# 1. Get significant genes from your DE results
res_df <- as.data.frame(res)
sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
sig_gene_ids <- rownames(sig_genes)
up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]

# significant genes
length(sig_gene_ids)   

# 2. Strip Ensembl version numbers (e.g. "ENSG00000063176.16" -> "ENSG00000063176")
gene_ids_clean <- sub("\\..*", "", sig_gene_ids)

# 3. Convert Ensembl IDs -> Entrez IDs (KEGG requires Entrez)
gene_map <- bitr(gene_ids_clean,
                 fromType = "ENSEMBL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

head(gene_map)
nrow(gene_map)          # how many successfully mapped (some will fail to map - that's normal)

# 4. KEGG pathway enrichment
kegg_result <- enrichKEGG(gene = gene_map$ENTREZID,
                          organism = "hsa",     # hsa = Homo sapiens
                          pvalueCutoff = 0.05)

# 5. View results as a dataframe
kegg_df <- as.data.frame(kegg_result)
head(kegg_df)

# 6. Visualize
dotplot(kegg_result, showCategory = 15)
barplot(kegg_result, showCategory = 15)