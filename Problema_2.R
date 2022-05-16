rm(list=ls(all=TRUE))
options(scipen = 100)
set.seed(2015790003)
n = 10000
p = 1000

H = matrix(nrow = n, ncol = n)
for(i in 1:n) {
       H[i,1]=1/sqrt(n) }

for (i in 1:n-1) {
       for (j in (i+1):n) {
        H[i,j]= 1/sqrt(j*(j-1))    
         }}

for (i in 2:n) {
  H[i,i]= -(i-1)/sqrt(i*(i-1))
  }

for (i in 1:n) {
  for (j in 1:n) {
    if(is.na( H[i,j]))
    {
    H[i,j]= 0    
    }}}
dim(H)

beta      <- rep(0,p) 
beta[1:5] <- 1:5

X= matrix(rnorm(H[,1:p]), nrow = n, ncol = p)
dim(X)
X <- scale(X)
Xb <- X%*%beta
e= rnorm(n)

Y <- X%*%beta+e
Y <- Y-mean(Y) 
plot(cor(X,Y),xlab="j (variables)",ylab="Coeff of Cor(Y,X_j)",main="Sample correlations",cex=2)

ols <- lm(Y~X)
beta_ols <- ols$coef[-1]
coef_ols=beta_ols[1:10]
summary(ols)
round(summary(ols)$coef[summary(ols)$coef[,4] < .06, 1:4], digits = 7)
##### LASSO #####
library(knitr)
library(lars)
lasso  <- lars(X,Y)
plot(lasso)

betas    <- lasso$beta
df       <- lasso$df # degree of freedom
MSE      <- lasso$RSS/n # mean squared error
bic      <- log(n)*df+n*log(MSE)

kable(cbind(df,MSE,bic))
bestb    <- which.min(bic)
bestb

plot(bic,cex=2)
points(bestb,bic[bestb],pch=19,cex=2, col="chocolate1")

bestMSE   <- which.min(MSE)
plot(MSE,cex=2)
points(bestMSE,MSE[bestMSE],pch=19,cex=2, col="red")

beta_lasso <- betas[bestb,]
kable(head(beta_lasso, n = 10), col.names = "Coefficients")
coef_lasso <- head(beta_lasso, n = 10)
##### CV #####
library(glmnet)
glmnet <- glmnet(X, Y, alpha = 1, lambda = lasso$lambda,standardize=FALSE) 
plot(glmnet)
cv <- cv.glmnet(X,Y, alpha = 1, lambda = lasso$lambda,standardize=FALSE)
best_lambda = cv$lambda.min
best_lambda
plot(cv)

glmnet_best_lambda=glmnet(X, Y, alpha = 1, lambda = best_lambda , standardize=FALSE)
beta_glmnet_best_lambda=coef(glmnet_best_lambda)[-1]
kable(head(beta_glmnet_best_lambda,n=10 ),col.names = "Coefficients")
coef_glmnet_best_lambda=head(beta_glmnet_best_lambda, n=10)

kable(cbind(coef_ols, coef_lasso, coef_glmnet_best_lambda))

#incremento la sd

e_sd2=rnorm(n, sd=2)
Y_sd2 <- X%*%beta+e_sd2
Y_sd2 <- Y_sd2-mean(Y_sd2) 

ols_sd2 <- lm(Y_sd2~X)
beta_ols_sd2= ols_sd2$coef[-1]
summary(ols_sd2)
round(summary(ols_sd2)$coef[summary(ols_sd2)$coef[,4] < .06, 1:4], digits = 10)
coef_ols_sd2=beta_ols_sd2[1:10]
lasso_sd2 <- lars(X, Y_sd2)
plot(lasso_sd2)

      
betas_sd2    <- lasso_sd2$beta
df_sd2       <- lasso_sd2$df # degree of freedom
MSE_sd2      <- lasso_sd2$RSS/n # mean squared error
bic_sd2      <- log(n)*df_sd2+n*log(MSE_sd2)

kable(cbind(df_sd2,MSE_sd2,bic_sd2))
bestb_sd2    <- which.min(bic_sd2)
bestb_sd2

plot(bic_sd2,cex=2)
points(bestb_sd2,bic_sd2[bestb_sd2],pch=19,cex=2, col="red")

beta_lasso_sd2 <- betas_sd2[bestb_sd2,]
coef_lasso_sd2 <- head(beta_lasso_sd2, n = 10)

glmnet_sd2 <- glmnet(X, Y_sd2, alpha = 1, lambda = lasso_sd2$lambda,standardize=FALSE) 
plot(glmnet_sd2)
cv_sd2 <- cv.glmnet(X,Y_sd2, alpha = 1, lambda = lasso_sd2$lambda,standardize=FALSE)
best_lambda_sd2 = cv_sd2$lambda.min
best_lambda_sd2
plot(cv_sd2)

glmnet_best_lambda_sd2=glmnet(X, Y, alpha = 1, lambda = best_lambda_sd2 , standardize=FALSE)
beta_glmnet_best_lambda_sd2=coef(glmnet_best_lambda_sd2)[-1]
coef_glmnet_best_lambda_sd2=head(beta_glmnet_best_lambda_sd2, n=10)

kable(cbind(coef_ols_sd2, coef_lasso_sd2, coef_glmnet_best_lambda_sd2))

kable(cbind(coef_ols,coef_ols_sd2))
kable(cbind(coef_lasso,coef_lasso_sd2))
kable(cbind(coef_glmnet_best_lambda, coef_glmnet_best_lambda_sd2))                               

MSE_ols <- mean((beta-beta_ols)^2)

MSE_lasso <- mean((beta-beta_lasso)^2)

MSE_glmnet <- mean((beta-beta_glmnet_best_lambda)^2)

MSE_ols2 <- mean((beta-beta_ols_sd2)^2)

MSE_lasso2 <- mean((beta-beta_lasso_sd2)^2)
MSE_glmnet2 <- mean((beta-beta_glmnet_best_lambda_sd2)^2)

kable(rbind(MSE_ols, MSE_lasso, MSE_glmnet, MSE_ols2, MSE_lasso2, MSE_glmnet2), col.names = "Mean Squared Error")

