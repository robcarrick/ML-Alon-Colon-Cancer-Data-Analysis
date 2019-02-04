### Assignment 4: Alon data set exploration ###
### Import Data ####
load('alon.rDa')
data <- cbind.data.frame(x,y)

### Feature correlation and p-value visualisations ####
# Corrleation histogram: histogram of feature correlations
res <- cor(x)
res[upper.tri(res, diag=TRUE)] <- 0
corr_vec <- res[res!=0]
hist(corr_vec, 75, xlab="Correlation", main="Correlation Histogram")

# Five figure summary of the correlations
corr_summary <- quantile(corr_vec)
names(corr_summary) <- c("Min","Q1","Median","Q3","Max")

### Pairs plot of top 4 covariates

# First run binary logistic regression to find the p values
p_vals <- integer(ncol(x)) # Initialise the p value vector

for (i in 1:ncol(x)){
  # Run logistic regression
  log_reg_mod <- glm(y ~ x[ ,i], family='binomial')
  p_vals[i] <- coef(summary(log_reg_mod))[,4][2]
}

# Order the p values (keep the indices)
p_vals_ordered <- sort(p_vals,index.return=TRUE)

top_4 <- p_vals_ordered$ix[1:4] # take the four covariates with the smallest p-values

# Plot
pairs(x[,top_4],main = 'Alon Dataset: Top Four Covariates', 
      lower.panel=NULL, 
      pch = c(22,24)[unclass(y)], 
      bg = c("red", "green")[unclass(y)])

# Histogram of the p values
hist(p_vals, 50,
     xlab = "P-Values",
     main = "Histogram of P-Values")

