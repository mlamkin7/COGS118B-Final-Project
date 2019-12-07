# number of genes 1000, note group prob must sum to 1 (only impacts number of cells for each cell line so in our case assume equal uniformity)

library(splatter)
library(scater)
library(stats)
library(anocva)
library(apcluster)
library(mclust)

sim.groups <- splatSimulateGroups(batchCells = 1000,
                                  group.prob = c(0.26, 0.12, 0.10, 0.37,
                                                 0.15),
                                  de.prob = c(0.05, 0.1, 0.08, 0.06,
                                              0.2), de.facLoc = 0.1,
                                  de.facScale = 0.4,
                                  seed = 1)
sim.groups <- normalize(sim.groups)
tcounts.matrix <- t(logcounts(sim.groups))
propergroups <- sim.groups@colData$Group

# t-sne plot
plotTSNE(sim.groups, colour_by = "Group")

# Hierarchical Clustering
d <- dist(tcounts.matrix)
hc <- hclust(d)
plot(hc)
hc.clusters <- rect.hclust(hc, k = 5, border = "red")

# run spectral clustering
spect.matrix <- expSimMat(tcounts.matrix)
spec.clusters <- spectralClustering(spect.matrix, 5)

# run K-means
k.clusters <- kmeans(tcounts.matrix, 5)
k.clusters <- k.clusters$cluster

# UNABLE TO RUN GMM DUE TO MEMORY ISSUES FOR PRE PCA

# Run PCA and use first two PCs for clustering of points again
counts.transformed <- plotPCA(sim.groups)$data

# Hierarchical Clustering
d <- dist(counts.transformed)
thc <- hclust(d)
plot(thc)
thc.clusters <- rect.hclust(thc, k = 5, border = "red")

# run spectral clustering
spect.matrix <- expSimMat(counts.transformed)
tspec.clusters <- spectralClustering(spect.matrix, 5)

# run K-means
tk.clusters <- kmeans(counts.transformed, centers=5)
tk.clusters <- tk.clusters$cluster

# Run GMM
tfit <- Mclust(counts.transformed, G=5)
tgmm.clusters <- tfit$classification

# Make GMM, Spectral, and K-means plots
# K-means plot
library(useful)
plot.kmeans(tk.clusters, data=counts.transformed, title="K-means on PCA transofrmed scRNA Expression Data")

# Spectral plot
plot(counts.transformed, col=tspec.clusters, xlab="PCA 1", ylab="PCA 2")

# GMM Plot
plot.Mclust(tfit, what=c("classification"), xlab="PCA 1", ylab="PCA 2")

# Evaluate clustering by outsourcing the labels to python
# write actual labels
write(propergroups, file="true_groups.csv", ncolumns=1)

# 7 total pre PCA: hierarchical, kmeans, spectral
#         post PCA: hierarchical, kmeans, spectral, GMM
write(k.clusters, file="/home/michael/COGS118B/Pre_PCA_kmeans.csv", ncolumns=1)
write(spec.clusters, file="/home/michael/COGS118B/Pre_PCA_spectral.csv", ncolumns=1)
write(tk.clusters, file="/home/michael/COGS118B/Post_PCA_kmeans.csv", ncolumns=1)
write(tspec.clusters, file="/home/michael/COGS118B/Post_PCA_spectral.csv", ncolumns=1)
write(tgmm.clusters, file="/home/michael/COGS118B/Post_PCA_GMM.csv", ncolumns=1)

fnlist <- function(x, fil){ z <- deparse(substitute(x))
cat(z, "\n", file=fil)
nams=names(x) 
for (i in seq_along(x) ){ cat(nams[i], "\t",  x[[i]], "\n", 
                              file=fil, append=TRUE) }
}
fnlist(hc.clusters, "/home/michael/COGS118B/Pre_PCA_hierarchical.csv")
fnlist(thc.clusters, "/home/michael/COGS118B/Pre_PCA_hierarchical.csv")
