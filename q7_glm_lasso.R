### Assignment 4: Q7 GLM Lasso ###
library(glmnet)
library(caret)
### Import Data ####
load('C:/Users/robca/Desktop/University/Year 4, 2018/Semester 2/STAT4401/Assignment 4/alon.rDa')
data <- cbind.data.frame(x,y)

### Part B: Tune the lambda parameter with cv ####

set.seed(1)
cv_lasso <- cv.glmnet(x, y, 
                      family="binomial", 
                      type.measure="class",
                      nfolds=5)

lambda_optimal <- cv_lasso$lambda.min
print(lambda_optimal)
plot(cv_lasso)
  
### Train optimal lasso glm on all the data and find error rates

lasso_mod_full <- glmnet(x, y, family="binomial", lambda=lambda_optimal)

# Make predictions for the entire data set
lasso_preds <- predict(lasso_mod_full, x, type='class')

# Extract the error rates from the confusion matrix
conf_mat_lasso <- confusionMatrix(as.factor(lasso_preds), y)

### Part C: CV estimates for error rates ####
# Overal error --------------------------------------------------0

set.seed(1)
train_control <- trainControl(method="cv", number=5)

model_lasso <- train(x=x, y=y, 
               trControl=train_control, 
               method="glmnet",
               metric="Accuracy", 
               tuneGrid=expand.grid(.alpha = 1,
                                    .lambda = lambda_optimal)
               )
cv_error_lasso <- 1 - model_lasso$results$Accuracy

# "Normal" class error ------------------------------------------0

normal_acc <- function(data, lev=NULL, model=NULL){
  out <- posPredValue(data[, "pred"], data[, "obs"], positive=lev[1])
  names(out) <- "Prec"
  out
}

set.seed(1)
train_control <- trainControl(summaryFunction = normal_acc, method="cv", number=5)

normal_cv_lasso <- train(x=x, y=y, 
               trControl=train_control, 
               method="glmnet",
               metric="Prec",
               tuneGrid=expand.grid(.alpha = 1,
                                    .lambda = lambda_optimal)
)

normal_cv_error_lasso <- 1 - normal_cv_lasso$results$Prec

# "Tumour" class error ------------------------------------------0

tumour_acc <- function(data, lev=NULL, model=NULL){
  out <- posPredValue(data[, "pred"], data[, "obs"], positive=lev[2])
  names(out) <- "Prec"
  out
}

set.seed(1)
train_control <- trainControl(summaryFunction = tumour_acc, method="cv", number=5)

tumour_cv_lasso <- train(x=x, y=y, 
                         trControl=train_control, 
                         method="glmnet",
                         metric="Prec",
                         tuneGrid=expand.grid(.alpha = 1,
                                              .lambda = lambda_optimal)
)

tumour_cv_error_lasso <- 1 - tumour_cv_lasso$results$Prec

### Part G: Plot of loss against lambda ####

### Plot of errors against lambda

set.seed(1)
cv_lasso <- cv.glmnet(x, y, 
                      family="binomial", 
                      type.measure="class",
                      nfolds=5)
# Plot the misclassification error against log(lambda)
plot(cv_lasso)

### Train the optimal lasso glm on all the data

optimal_lasso_glm <- glmnet(x, y,
                            family="binomial",
                            lambda = lambda_optimal, 
                            )

# Extract the non-zero coefficients
lasso_coef <- coef(optimal_lasso_glm, s = lambda_optimal)
lasso_covariates <- which(lasso_coef!=0) - 1 # -1 because of the intercjept term
lasso_coef_nonzero <- lasso_coef[which(lasso_coef!=0)]
print(lasso_coef_nonzero)

# Format as a dataframe
coefs <- as.data.frame(matrix(0, nrow=length(lasso_coef_nonzero), ncol=2))
colnames(coefs) <- c("Gene Number","Est Coefficient")
coefs[ ,1] <- rownames(lasso_coef)[which(lasso_coef!=0)]
coefs[ ,2] <- lasso_coef_nonzero
