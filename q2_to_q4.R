### Assignment 4: Alon data set exploration ###
### Import Data ####
load('C:/Users/robca/Desktop/University/Year 4, 2018/Semester 2/STAT4401/Assignment 4/alon.rDa')
data <- cbind.data.frame(x,y)

### Question 2: Visualisation ####

hist(x[ ,200])
plot(x[ ,1], x[ ,2])
pairs(x[ ,1:6], cex=0.5, pch=20)

plot(x[ ,625], as.numeric(y)-1)
hist(as.numeric(lapply(x, min)), 15)
hist(as.numeric(lapply(x, max)), 15)

# Corrleation histogram
res <- cor(x)
res[upper.tri(res, diag=TRUE)] <- 0
corr_vec <- res[res!=0]
hist(corr_vec, 75, xlab="Correlation", main="Correlation Histogram")

# Five figure summary of the correlations
corr_summary <- quantile(corr_vec)
names(corr_summary) <- c("Min","Q1","Median","Q3","Max")

### Pairs plot of top 4 covariates

# First run logisitc regression to find the p values
p_vals <- integer(ncol(x)) # Initialise the p value vector

for (i in 1:ncol(x)){
  # Run logistic regression
  log_reg_mod <- glm(y ~ x[ ,i], family='binomial')
  p_vals[i] <- coef(summary(log_reg_mod))[,4][2]
}

# Order the p values (keep the indices)
p_vals_ordered <- sort(p_vals,index.return=TRUE)

top_4 <- p_vals_ordered$ix[1:4]

# Plot
pairs(x[,top_4],main = 'Alon Dataset: Top Four Covariates', 
      lower.panel=NULL, 
      pch = c(22,24)[unclass(y)], 
      bg = c("red", "green")[unclass(y)])

# Histogram of the p values
hist(p_vals, 50,
     xlab = "P-Values",
     main = "Histogram of P-Values")

### Question 3: PCA (With FactoMineR) ####
library(FactoMineR)

#alon_pca <- PCA(X = x, ncp = 61, graph=T)

alon_pca <- PCA(X = data, ncp = 61, scale.unit = FALSE, 
                quali.sup = 2001, graph=T)

plot.PCA(alon_pca, axes=c(1,2), choix="ind", habillage=2001,
         legend=list(x="topright"))

# Plot the % of variance wrt principal component
plot(alon_pca$eig[ ,2], 
     pch=1,
     cex=0.5,
     xlab='Principal Component', 
     ylab='% of Variance Explained',
     main='PCA: % of Variance Explained')
lines(alon_pca$eig[ ,2], lty='dotted')

# Plot the cumulative % of variance wrt principal component
plot(alon_pca$eig[ ,3], 
     pch=1,
     cex=0.5,
     xlab='Principal Component', 
     ylab='Cumulative % of Variance Explained',
     main='PCA: Cumulative % of Variance Explained')
lines(alon_pca$eig[ ,3], lty='dotted')

### Store the relevant outputs
# 
# pca_eig <- alon_pca$eig # Eigenvalues for each PC and the indiv and cum % of var explained
# pca_PCs <- alon_pca$var$coord # All 61 principal components 
# 
# # Directions of the principal components
# alon_pca_directions <- matrix(0, nrow=nrow(alon_pca$var$coord),
#                               ncol=ncol(alon_pca$var$coord))
# rownames(alon_pca_directions) <- rownames(alon_pca$var$coord)
# colnames(alon_pca_directions) <- colnames(alon_pca$var$coord)
# 
# for (i in 1:ncol(alon_pca_directions)){
#   alon_pca_directions[ ,i] <- alon_pca$var$coord[ ,i] /
#                               sqrt(sum(alon_pca$var$coord[ ,i]^2))
# }
# 
# pca_directions <- alon_pca_directions
# pca_x_map <- alon_pca$ind$coord # Observations of x mapped to principal components
# 
# save(pca_eig, pca_PCs, pca_directions, pca_x_map, 
#      file = "C:/Users/robca/Desktop/University/Year 4, 2018/Semester 2/STAT4401/Assignment 4/pca_output.rda")
# 
# load('C:/Users/robca/Desktop/University/Year 4, 2018/Semester 2/STAT4401/Assignment 4/pca_output.rda')


### Question 4: Single Variable Analysis ####

# Use the ordered p-value list from question 2

# Find the value of k
fdr <- 0.02

ind <- integer()
j <- 1
for (i in 1:ncol(x)){
  if (p_vals_ordered$x[i] <= ((i/ncol(x))*fdr)){
    ind[j] <- i
    j <- j + 1
  }
}
k <- max(ind)

# Now take the k covariates with the smallest p-values
sig_covariates <- p_vals_ordered$ix[1:k]
vals <- p_vals_ordered$x[1:k]

# Get Adjusted P Values
p_vals_adj <- integer(ncol(x))
for (i in 1:ncol(x)){
  p_vals_adj[i] <- p_vals_ordered$x[i] * ncol(x) / i
}

### Plot

# Find the threshold line
threshold <- integer(ncol(x))
for (i in 1:ncol(x)){
  threshold[i] <- fdr * i/ncol(x)
}

plot(p_vals_ordered$x, type="l", col='blue',
     xlab='Order of P Value', ylab='P Value',
     main='Benjamini-Hochberg FDR Control')
lines(threshold, type="l", col="red")
#legend(5,0.0025,legend = c('P Value','Threshold'),col = c("blue","red"),
#       lty=c(1,1))
legend(5,0.8,legend = c('P Value','Threshold'),col = c("blue","red"),
       lty=c(1,1))

### List of significant genes with associated threshold values

sig_genes <- as.data.frame(matrix(0,nrow=length(sig_covariates),ncol=3))
colnames(sig_genes) <- c("Gene Number","P Value","FDR Threshold")
sig_genes[ ,1] <- sig_covariates
sig_genes[ ,2] <- p_vals_ordered$x[1:length(sig_covariates)]
sig_genes[ ,3] <- threshold[1:length(sig_covariates)]

## Save as an rda file
# save(sig_genes,
#      file = "C:/Users/robca/Desktop/University/Year 4, 2018/Semester 2/STAT4401/Assignment 4/sig_genes.rda")
load("C:/Users/robca/Desktop/University/Year 4, 2018/Semester 2/STAT4401/Assignment 4/sig_genes.rda")
