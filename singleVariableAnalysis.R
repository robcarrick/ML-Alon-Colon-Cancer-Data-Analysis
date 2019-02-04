### Import Data ####
load('alon.rDa')
data <- cbind.data.frame(x,y)

### Single variable analysis: false discovery rate (FDR) is controlled ####

# Run binary logistic regression with each covariate separately (2000 models) and extract p-values
p_vals <- integer(ncol(x)) # Initialise the p value vector
for (i in 1:ncol(x)){
  # Run logistic regression
  log_reg_mod <- glm(y ~ x[ ,i], family='binomial')
  p_vals[i] <- coef(summary(log_reg_mod))[,4][2]
}
# Order the p values (keep the indices)
p_vals_ordered <- sort(p_vals,index.return=TRUE)

### Implement Benjamini-Hochberg correction (controlling the FDR)
# Find the value of k
fdr <- 0.02 # control the FDR to be at most 2%
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

# Plot gene order by p-value and include the FDR control threshold line
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

### List of the k significant genes with associated p-values and threshold values
sig_genes <- as.data.frame(matrix(0,nrow=length(sig_covariates),ncol=3))
colnames(sig_genes) <- c("Gene Number","P Value","FDR Threshold")
sig_genes[ ,1] <- sig_covariates
sig_genes[ ,2] <- p_vals_ordered$x[1:length(sig_covariates)]
sig_genes[ ,3] <- threshold[1:length(sig_covariates)]
