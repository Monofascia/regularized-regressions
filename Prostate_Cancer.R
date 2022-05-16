rm(list = ls(all=TRUE))
library(ISLR)
library(glmnet)
library(dplyr)
library(tidyr)
library(readr)
library(knitr)
library(psych)

options(scipen = 10)
prostate <- read_delim("C:/Users/Patrizio/Desktop/TUTTO/Ud'A/CLEBA/statistical learning/HOMEWORK REGULARIZED REGRESSION/prostate.data", 
                       "\t", escape_double = FALSE, col_types = cols(X1 = col_skip()), 
                       trim_ws = TRUE)
prostate$svi= as.logical(prostate$svi)

variabili <- c("lpsa","lcavol", "lweight", "age", "lbph", "svi", "lcp",
               "gleason", "pgg45")
variabili
descrizione<- c("level of prostate-specific antigen", 
                "log(cancer volume)", "log(prostate weight)", 
                "age", "log(benign prostatic hyperplasia amount)", 
                "seminal vesicle invasion", "log(capsular penetration)", 
                "Gleason score", "percentage Gleason scores 4 or 5")
descrizione
kable(cbind(variabili, descrizione))

kable(head(prostate))

# train
train <- subset(prostate, train == "TRUE")
train_FREE <- train[,-10] #elimina colonna train
x_train_matrix <- train_FREE %>% select(-lpsa) %>% as.matrix() 
y_train <-  train_FREE %>% 
  select(lpsa) %>% 
  scale(center = TRUE, scale = FALSE) %>% 
  as.matrix() 

#test
test <- subset(prostate, train == "FALSE")
test_FREE <- test[,-10] #elimino colonna train
x_test_matrix <- test_FREE %>% select(-lpsa) %>% as.matrix() 
y_test <- test_FREE %>% 
  select(lpsa) %>% 
  scale(center = TRUE, scale = FALSE) %>% 
  as.matrix() 

########### analisi esplorativa #################
par(mfrow = c(1,2))
plot( train_FREE$lcavol , 
      train_FREE$lpsa , 
      cex = 0.5 , pch=19  , 
      main = "SCATTERPLOT Y ed lcavol " , xlab = "X=lcavol" , ylab = "Y=lpsa")
plot( train_FREE$lweight , 
      train_FREE$lpsa ,
      cex = 0.5 , 
      pch=19  , 
      main = "SCATTERPLOT Y ed lweight " , xlab = "X=lweight" , ylab = "Y=lpsa")
par(mfrow = c(1,2))
plot( train_FREE$age , 
      train_FREE$lpsa ,
      cex = 0.5 , 
      pch=19  , 
      main = "SCATTERPLOT Y e age" , xlab = "X=age" , ylab = "Y=lpsa")
plot( train_FREE$svi , 
      train_FREE$lpsa ,
      cex = 0.5 , 
      pch=19  , 
      main = "SCATTERPLOT Y e svi " , xlab = "X=svi" , ylab = "Y=lpsa")


########### OLS ###########
 
lm_TRAIN<- lm(lpsa~., data = train_FREE)
summary(lm_TRAIN)

summary.lm(lm_TRAIN, signif.stars = TRUE)


MSE_lm_TRAIN <- mean((train_FREE$lpsa - predict(lm_TRAIN))^2)
prediction_of_lm_TRAIN_on_TEST <- predict(lm_TRAIN, newdata = test_FREE)
MSE_lm_TRAIN_on_TEST <- mean((test_FREE$lpsa - prediction_of_lm_TRAIN_on_TEST)^2)
kable(cbind(MSE_lm_TRAIN_on_TEST, MSE_lm_TRAIN), 
      col.names = c("MSE su TEST","MSE su TRAIN"))
library(car)
kable(vif(lm_TRAIN), col.names = "Mean squared error")
vif_value <- vif(lm_TRAIN)
barplot(vif_value, main = "VIF Values",
        col = rainbow(8), ylim = c(0,6), space = 0.5)
abline(v = 5, lwd = 3, lty = 2)

# b
lm_lcavol <- lm(lpsa~lcavol, data =train_FREE) #lcavol
summary(lm_lcavol)
lm_lcavol_lweight<- lm(lpsa~lcavol+lweight, data = train_FREE) #+lweight
summary(lm_lcavol_lweight)
lm_lcavol_lweight_svi<- lm(lpsa~lcavol+lweight+svi, data = train_FREE) #+svi
summary(lm_lcavol_lweight_svi)
lm_lcavol_lweight_svi_lbph<- lm(lpsa~lcavol+lweight+svi+lbph, data = train_FREE) #+lbph
summary(lm_lcavol_lweight_svi_lbph)
lm_lcavol_lweight_svi_lbph_lcp<- lm(lpsa~lcavol+lweight+svi+lbph+lcp, data = train_FREE) #+lcp
summary(lm_lcavol_lweight_svi_lbph_lcp)
lm_lcavol_lweight_svi_lbph_lcp_pgg45<- lm(lpsa~lcavol+lweight+svi+lbph+lcp+pgg45, 
                                          data = train_FREE) #+pgg45
summary(lm_lcavol_lweight_svi_lbph_lcp_pgg45)
lm_lcavol_lweight_svi_lbph_lcp_pgg45_age<- lm(lpsa~lcavol+lweight+svi+lbph+lcp+pgg45+age,
                                              data = train_FREE) #+age
summary(lm_lcavol_lweight_svi_lbph_lcp_pgg45_age)

# BIC
BIC_OLS<-BIC(lm_TRAIN, lm_lcavol,lm_lcavol_lweight,lm_lcavol_lweight_svi, 
    lm_lcavol_lweight_svi_lbph, lm_lcavol_lweight_svi_lbph_lcp, 
    lm_lcavol_lweight_svi_lbph_lcp_pgg45, lm_lcavol_lweight_svi_lbph_lcp_pgg45_age )
kable(BIC_OLS)
min(BIC_OLS)

kable(coef(lm_lcavol_lweight), col.names = "Linear model")

# MSE su Test e Train
MSE_best_lm <- mean((train_FREE$lpsa - predict(lm_lcavol_lweight))^2)
MSE_best_lm_on_TEST<- mean((test_FREE$lpsa - predict(lm_lcavol_lweight, newdata = test_FREE))^2)
kable(cbind(MSE_best_lm_on_TEST, MSE_best_lm))

kable(cbind(MSE_best_lm_on_TEST, MSE_best_lm), 
      col.names = c("MSE su Test", "MSE su Train"))

R2_TRAIN<-summary(lm_TRAIN)$r.squared
R2_best<-summary(lm_lcavol_lweight)$r.squared
kable(cbind(R2_TRAIN,R2_best))

############ RIDGE ############
set.seed(1)
grid = 10^seq(10, -2, length = 100)
ridge_mod <- glmnet(x_train_matrix , y_train, alpha = 0, lambda = grid) 
dim(coef(ridge_mod))
plot(ridge_mod)
cv.out = cv.glmnet(x_train_matrix, y_train, alpha = 0, lambda = grid) 
best_lambda = cv.out$lambda.min
best_lambda
plot(cv.out)

ridge_mod_TOP = glmnet(x_train_matrix , y_train, alpha = 0, lambda = best_lambda)
MSE_Ridge_cv <- mean((y_train - predict(ridge_mod_TOP, newx = x_train_matrix))^2)
MSE_Ridge_cv_TEST <- mean((y_test - predict(ridge_mod_TOP, newx = x_test_matrix))^2)

############# BIC for best lambda in Ridge #######################


n_Ridge <- ridge_mod$nobs
n_Ridge
MSE_Ridge <- deviance(ridge_mod)/n_Ridge
MSE_Ridge
df_Ridge <- ridge_mod$df
df_Ridge
BIC <- log(n_Ridge)*df_Ridge + n_Ridge*log(MSE_Ridge) 
BIC
min(BIC)
which.min(BIC)
ridge_mod$lambda[100]

ridge_mod_BIC <- glmnet(x_train_matrix, y_train, alpha = 0, lambda = grid[100])
MSE_Ridge_bic <- mean((y_train - predict(ridge_mod_BIC, newx = x_train_matrix))^2)
MSE_Ridge_bic_TEST <- mean((y_test - predict(ridge_mod_BIC, newx = x_test_matrix))^2)

# MSE su Test e Train

kable(cbind(MSE_Ridge_bic, MSE_Ridge_cv))
kable(cbind(MSE_Ridge_bic_TEST, MSE_Ridge_cv_TEST))

############ LASSO ###########


lasso_mod = glmnet(x_train_matrix, y_train, alpha = 1,lambda = grid)
plot(lasso_mod)
cv.out.lasso = cv.glmnet(x_train_matrix, y_train, alpha = 1, lambda = grid)
best_lambda_lasso = cv.out.lasso$lambda.min
best_lambda_lasso
plot(cv.out.lasso)

lasso_mod_TOP = glmnet(x_train_matrix, y_train, alpha = 1,lambda = best_lambda_lasso)
MSE_Lasso_cv <- mean((y_train - predict(lasso_mod_TOP, newx = x_train_matrix))^2)
MSE_Lasso_cv_TEST <- mean((y_test - predict(lasso_mod_TOP, newx = x_test_matrix))^2)


####### BIC for best lambda in LASSO ########## 

n_lasso <- lasso_mod$nobs
MSE_lasso <- deviance(lasso_mod)/n_lasso
df_lasso <- lasso_mod$df
BIC_lasso <- log(n_lasso)*df_lasso + n_lasso*log(MSE_lasso) 
min(BIC_lasso)
which.min(BIC_lasso)
bic_top <- which.min(BIC_lasso)
lambda_lasso_bic <- lasso_mod$lambda[93]
lambda_lasso_bic
plot(BIC_lasso,xlab = "Modello stimato", cex=2) 
points(bic_top,BIC_lasso[bic_top],pch=19,cex=2)  



###

lasso_mod_BIC = glmnet(x_train_matrix, y_train, alpha = 1, lambda =  lambda_lasso_bic )
MSE_Lasso_bic <- mean((y_train - predict(lasso_mod_BIC, newx = x_train_matrix))^2)
MSE_Lasso_bic_TEST <- mean((y_test - predict(lasso_mod_BIC, newx = x_test_matrix))^2)

# MSE su Test e Train
kable(cbind(MSE_Lasso_bic, MSE_Lasso_cv))
kable(cbind(MSE_Lasso_bic_TEST, MSE_Lasso_cv_TEST))


# MSE finale solo su TEST con tutti i modelli
kable(rbind(MSE_lm_TRAIN_on_TEST, 
            MSE_best_lm_on_TEST, 
            MSE_Ridge_cv_TEST, 
            MSE_Lasso_bic_TEST), 
      col.names = "Mean Squared Error")

# confronto coefficienti

beta_lm_TRAIN <- coef(lm_TRAIN)
beta_lm_BIC <- coef(lm_lcavol_lweight)
beta_ridge_CV <- coef.glmnet(ridge_mod_TOP) # quello con cv
beta_lasso_BIC <- coef.glmnet(lasso_mod_BIC)


matriceRis <- matrix(,nrow=9,ncol=4, 
                     dimnames = list(c("Intercept","lcavol", "lweight", "age", "lbph", "svi", "lcp", "gleason", "pgg45"),
                                     c("beta lm TRAIN" , "beta lm BIC","beta ridge CV","beta lasso BIC")
                                     )) 
                                     
                     


for (i in 1:9) {
  
  matriceRis[i,1] <- beta_lm_TRAIN[i]
  
}

for (i in 1:3) {
  
  matriceRis[i,2] <- beta_lm_BIC[i]
  
  
}

for (i in 1:9) {
  
  matriceRis[i,3] <- beta_ridge_CV[i]
  
}

for (i in 1:9) {
  
  matriceRis[i,4] <- beta_lasso_BIC[i]
  
}

for (i in 4:9) {
  
 matriceRis[i,2] <- 0.000000000
  
}


kable(matriceRis)


############## FINISH ################
  