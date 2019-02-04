### Import Data ####
load('alon.rDa')
data <- cbind.data.frame(x,y)

### PCA (With FactoMineR) ####
library(FactoMineR)
alon_pca <- PCA(X = data, ncp = 61, scale.unit = FALSE, 
                quali.sup = 2001, graph=T)
plot.PCA(alon_pca, axes=c(1,2), choix="ind", habillage=2001,
         legend=list(x="topright"))

### Extract relevant outputs
pca_eig <- alon_pca$eig # Eigenvalues for each PC and the indiv and cum % of var explained
pca_PCs <- alon_pca$var$coord # All 61 principal components 
pca_x_map <- alon_pca$ind$coord # Observations in x mapped to the linear sub-space spanned by the principal components

# Plot the % of variance wrt principal component
plot(pca_eig[ ,2], 
     pch=1,
     cex=0.5,
     xlab='Principal Component', 
     ylab='% of Variance Explained',
     main='PCA: % of Variance Explained')
lines(pca_eig[ ,2], lty='dotted')

# Plot the cumulative % of variance wrt principal component
plot(pca_eig[ ,3], 
     pch=1,
     cex=0.5,
     xlab='Principal Component', 
     ylab='Cumulative % of Variance Explained',
     main='PCA: Cumulative % of Variance Explained')
lines(pca_eig[ ,3], lty='dotted')

