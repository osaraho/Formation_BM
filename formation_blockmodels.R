rm(list=ls())
#---------------------------------------------
### datasets
#---------------------------------------------

load('fungi_tree_data.Rdata')
ls()

#---------------------------------------------
#### useful packages 
#---------------------------------------------
# install.packages("blockmodels")
library(blockmodels)
library(igraph)
library(alluvial)
library(ggplot2)
source('function_for_blockmodels.R')

#---------------------------------------------
#### SBM for a binary adjacency matrix
#-------------------------------------------

# data
plotMatrix(Mat = tree_bin,rowFG = 'tree', colFG  = 'tree')

# SBM for bernoulli 
if (!exists("sbm.tree_bin")) {
  sbm.tree_bin <- BM_bernoulli("SBM_sym",tree_bin)
  sbm.tree_bin$estimate()
}


#save(sbm.tree_bin,file='res_sbm_tree_bin.Rdata')
#load(file='res_sbm_tree_bin.Rdata')

# optimal number of blocks / cluster
Q = which.max(sbm.tree_bin$ICL) 
Q

#- extract the estimated parameters
paramEstimSBM <- extractParamBM(sbm.tree_bin,Q)
paramEstimSBM$pi
paramEstimSBM$alpha
paramEstimSBM$Z

#- rearrange the matrix followwing blocks
plotMatrix(tree_bin,'tree','tree', fileNameSave = NULL, clustering = list(row = paramEstimSBM$Z))

#- macropscopic view of the network
G <- graph_from_adjacency_matrix(paramEstimSBM$alpha, mode = c("undirected"), weighted = TRUE, diag = TRUE)
plot.igraph(G,vertex.size = paramEstimSBM$pi*100,edge.width= abs(E(G)$weight)*3,vertex.color = 1:Q, layout = layout_nicely)

#-  composition of the clusters
lapply(1:Q,function(q){tree_list[paramEstimSBM$Z == q]})

#---------------------------------------------
#### Poisson SBM for a weighted adjacency matrix
#-------------------------------------------

# data
plotMatrix(Mat = tree,rowFG = 'tree', colFG  = 'tree')


# SBM for Poisson
if (!exists("sbm.tree")) {
  sbm.tree <- BM_poisson("SBM_sym",tree)
  sbm.tree$estimate()
}
#save(sbm.tree,file='res_sbm_tree.Rdata')
#load(file='res_sbm_tree.Rdata')


# optimal number of blocks / cluster
Q = which.max(sbm.tree$ICL)
Q

# Extract parameters 
paramEstimSBMPoisson <- extractParamBM(sbm.tree,Q)
paramEstimSBMPoisson$pi
paramEstimSBMPoisson$lambda
paramEstimSBMPoisson$Z

# rearrange the matrix followwing blocks
plotMatrix(tree,'tree','tree', fileNameSave = NULL, clustering = list(row = paramEstimSBMPoisson$Z))

# macropscopic view of the network
G <- graph_from_adjacency_matrix(paramEstimSBMPoisson$lambda, mode = c("undirected"), weighted = TRUE, diag = TRUE)
plot.igraph(G,vertex.size = paramEstimSBMPoisson$pi * 200,edge.width = abs(E(G)$weight) / 2,vertex.color = 1:Q,layout = layout_nicely)
lapply(1:Q,function(q){tree_list[paramEstimSBMPoisson$Z == q]})


# Compare clustering with alluvial plots
A <- as.data.frame(table(paramEstimSBM$Z,paramEstimSBMPoisson$Z))
colnames(A) = c('SBM Bern',"SBM Poisson","Freq")
w   <- which(A$Freq != 0)
A <- A[w,]
alluvial(A[,c(1,2)],freq = A$Freq)


#---------------------------------------------
#### Poisson SBM with covariates
#---------------------------------------------


# SBM for Poisson
if (!exists("sbm.cov")) {
  sbm.cov <- BM_poisson_covariates("SBM",tree, ListVar)
  sbm.cov$estimate()
}

#save(sbm.cov,file='res_sbm_cov.Rdata')
#load(file='res_sbm_cov.Rdata')


# optimal number of blocks / cluster
Q = which.max(sbm.cov$ICL)
Q

# Extract parameters
paramEstimSBMPoissonCov <- extractParamBM(sbm.cov,Q)
paramEstimSBMPoissonCov$pi
paramEstimSBMPoissonCov$lambda
paramEstimSBMPoissonCov$Z
paramEstimSBMPoissonCov$beta  # effect of the covariates

# compare clusterings
B <- as.data.frame(table(paramEstimSBM$Z,paramEstimSBMPoisson$Z,paramEstimSBMPoissonCov$Z))
colnames(B) = c('SBM Bern',"SBM Poisson","SBM Cov","Freq")
w   <- which(B$Freq != 0)
B <- B[w,]
alluvial(B[,c(1,2,3)],freq = B$Freq)

# Macroscopic view 
G <- graph_from_adjacency_matrix(paramEstimSBMPoissonCov$lambda, mode = c("undirected"), weighted = TRUE, diag = TRUE)
plot.igraph(G,vertex.size=paramEstimSBMPoissonCov$pi*100,edge.width=sqrt(abs(E(G)$weight)),vertex.color = 1:Q, layout = layout_nicely)



#---------------------------------------------
####  Latent block model 
#---------------------------------------------

if (!exists("lbm")) {
  lbm <- BM_bernoulli("LBM",as.matrix(fungi_tree))
  lbm$estimate()
}

# extract better model index
Q = which.max(lbm$ICL)
Q
# extract optimal number of clusters
paramEstimLBM <- extractParamBM(lbm,Q)
paramEstimLBM$Q

# extract parameters
paramEstimLBM$piRow
paramEstimLBM$piCol
paramEstimLBM$alpha
paramEstimLBM$ZRow
paramEstimLBM$ZCol


# rearrange the matrix followwing blocks
plotMatrix(fungi_tree,'fungi','tree', fileNameSave = NULL, clustering = list(row = paramEstimLBM$ZRow,col = paramEstimLBM$ZCol))


# macropscopic view
G <- graph_from_incidence_matrix(paramEstimLBM$alpha, weighted = TRUE)
plot(G,vertex.size = c(paramEstimLBM$piRow * 100, paramEstimLBM$piCol * 100), vertex.shape = c("circle", "square")[V(G)$type + 1],        edge.width=abs(E(G)$weight*2),vertex.color=1:Q, layout=layout.bipartite)


