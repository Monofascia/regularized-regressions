#####################################################################################################
####################### HOMEWORK REGRESSIONE RIDGE E LASSO: PROBLEMA 3 ##############################
#####################################################################################################

# PULIZIA AMBIENTE DI LAVORO

rm(list=ls(all=TRUE))

#carico le librerie

library(glmnet) 
library(lars)

# IMPOSTO IL SEED

set.seed(2015790003)

# CREAZIONE DELLA MATRICE QUADRATA  H (n x n)

n = 6
p = 3

H = matrix(nrow = n, ncol = n)

# ciclo for 1: setta i valori della prima colonna della matrice H ad 1/sqrt(n)

for(i in 1:n) {
  
  H[i,1]=1/sqrt(n) 
  
  }

# ciclo for 2: setta i valori al di sopra della diagonale principale a 1/sqrt( n(n-1) )

for (i in 1:n-1) {
  
  for (j in (i+1):n) {
    H[i,j]= 1/sqrt(j*(j-1))    
  }
}

# ciclo for 3: setta i valori della diagonale principale ad  -(n-1)/sqrt( n(n-1) )

for (i in 2:n) {
  
  H[i,i]= -(i-1)/sqrt(i*(i-1))

  }
# ciclo for 4: setta i restanti elementi al di sotto della diagonale principale a zero 

for (i in 1:n) {
  
  for (j in 1:n) {
    
    if(is.na( H[i,j])) {
      H[i,j]= 0    
    }
  }
}

# DIAMO UN OCCHIATA AD H...

#head.matrix(H, n = 6)[1:6,1:6] # mostra a video le prime 6 righe e colonne di H

# generiamo il vero vettore dei coefficienti di regressione
#come coefficenti ho preso solo i primi 2. Secondo me il prof non aveva modificato il testo
#perchè all'inizio doveva essere quadrata anche X e quindi 6*6 e beta sarebbe stato: 1 2 3 4 5 0
beta      <- rep(0,p) # ripeti 0 per n volte
beta[1:2] <- 1:2     # solo i primi 2 coef. sono diversi da zero

# MATRICE DI SIMULAZIONE X3

#n1 <-14
#p1 <- 7

X1= matrix(rnorm(H[,1:p]), nrow = n, ncol = p)  # prende le prime p colonne di H e genera una mart. X (nxp)!!!????
X1 <- scale(X1)

# SIMULIAMO  Y 

#beta1 <- beta[1:3]

Xb1 <- X1%*%beta
e1= rnorm(n)         # reridui
Y1 <- X1%*%beta+e1   # perchè rifà X%*%beta?
Y1 <- Y1-mean(Y1) 

# RAPPRESENTAZIONE GRAFICA DELLE CORRELAZIONI

plot(cor(X1,Y1),xlab="j (variables)",ylab="Coeff of Cor(Y,X_j)",main="Sample correlations",cex=2)

# regressione OLS

regOLS <- lm(Y1~X1)
summary(regOLS)
coefOLS <- regOLS$coefficients[-1]
coefOLS

cor(X1)

# ALGORITMO DI KILLING

# 1) INIZIAZIZZAZIONE DEL VETTORE DEI COEFFICIENTI STIMATI MEDIANTE UNA RIDGE REGRESSION

lambdaGrid <-10^seq(10, -2, length = 100) # griglia di valori di lambda

ridgeReg <- glmnet(X1, Y1, alpha = 0, lambda = lambdaGrid,) # regressione ridge

plot(ridgeReg , main = list( "Ridge Regression Coefficient Paths" , cex = 0.8 ) , xlab =expression(lambda) 
     , ylab = "Coefficients")

# procedura di cross-validazione per inizializzare il vettore dei coefficienti beta

cv_ridgeReg <- cv.glmnet(X1,Y1, alpha = 0 , lambda = lambdaGrid, grouped = FALSE )
best_lambda <- cv_ridgeReg$lambda.min
best_lambda

plot(cv_ridgeReg , main = list( "RIDGE Regression cross-validation error Paths" , cex = 0.8 ) )
log(best_lambda) # punto in corrispondenza della prima linea trattegiata

# inizializzo il vettore dei coefficienti mediante una regressione ridge ed utilizzando il miglior
# lambda ottenuto mediante procedura di cross-validazione

ridgeReg_TOP <- glmnet(X1 , Y1 , alpha = 0 , lambda = best_lambda)
coef_TOP <- ridgeReg_TOP$beta #vettore dei coefficienti ridge
coef_TOP

cbind(coefOLS , coef_TOP)

# PROCEDURA DI KILLING

# riscrivo coef_TOP come vettore

#coef_Killing <- array(coef_TOP)   # inizializzo il vettore dei coefficienti beta 

coef_Killing_t <- coef_TOP  # vettore dei coefficienti stimati durante l'iterazione corrente
coef_Killing_t=as.vector(coef_Killing_t)
coef_iniziali=coef_Killing_t # si mantiene il vettore di coefficienti iniziale 
lambdaK <- 3.1424402 #con il valore di bestlambda della cross validation non funziona.

valueConv <- 0.0000000000000000001

fermati <- TRUE
count=0
coef_Killing_tp <- coef_Killing_t 
  
#iterazioni<- matrix(ncol = 4, nrow = 6)
#colnames(iterazioni)<-c("n°iterazione", "beta1","beta2","beta3")
label0 <- "vettore dei coefficienti Ridge iniziale: "   
label1 <- "iterazione: "
label2 <- "coefficienti stimati: "
  
  
  
  
while( fermati==TRUE ) {
  
  for (k in 1: length(coef_Killing_t)) {
  
    Xk <- X1[,k]                                # vettore del k-esimo regressore
    Xk1 <- cbind(Xk)                            # incolonno i valori di Xk => lo uso per il calcolo di Ck
    XkRestanti <- X1[,-k]                       # matrice dei restanti k-1 regressori
    betaRestanti <- coef_Killing_t[-k]          # vettore dei restanti k-1 parametri
  
    ak <- 2*(sum((Xk)^2))                       # calcolo di ak
  
    residui <- Y1 - (XkRestanti%*%betaRestanti) # faccio Y - X*beta (per i rimanenti k-1 regressori/parametri)
    Ck <- 2*(sum(residui * Xk1)) 
    
    if(Ck < -lambdaK){                          # confronto Ck e lambda per stimare betaK                       
    
      coef_Killing_t[k] <- (Ck + lambdaK)/ak
  
    } else { 
    
      if(Ck > lambdaK){
    
        coef_Killing_t[k] <- (Ck - lambdaK)/ak
  
      } else {
    
        coef_Killing_t[k] <- 0
        
        fermati=FALSE
        break
  
      }
  
   }

  }
  
  sumTp <- ( sum( abs(coef_Killing_tp ) ) ) # somma dei valori assoluti dei coefficienti stimati nell'iterazione precedente
  sumT <- ( sum( abs(coef_Killing_t ) ) )   # somma dei valori assoluti dei coefficienti stimati nell'iterazione corrente
  
  if(count==0){
    print("*******************************************************************************************" , quote = FALSE)
    print(" " , quote = FALSE)
    print( label0 , quote = FALSE ) # stampo il vettore dei coefficienti di partenza (Ridge)
    print(" " , quote = FALSE)
    print( c(" " , coef_iniziali ) , quote = FALSE ) # stampo il vettore dei coefficienti di partenza (Ridge)
    print(" " , quote =FALSE)
    
  }
  
  if( abs(sumT-sumTp) <= valueConv ){       # verifico la condizione di convergenza
         
      fermati <- FALSE
      count <- count+1
      print( c(label1 , count) , quote = FALSE )
      print(" " , quote =FALSE)
      print( c(label2 , coef_Killing_t) , quote = FALSE ) # stampa i coefficienti dell'iterazione corrente
      print(" " , quote =FALSE)
      print("*******************************************************************************************" , quote = FALSE)
      print(" " , quote = FALSE)
     
    }else{
              
      fermati <- TRUE
      coef_Killing_tp <- coef_Killing_t   
      count=count+1
      print( c(label1 , count) , quote = FALSE)
      print(" " , quote = FALSE)
      print( c(label2 , coef_Killing_tp ) , quote = FALSE)
      print(" " , quote = FALSE)
      print("*******************************************************************************************" , quote = FALSE)
      print(" " , quote = FALSE)
  }
  
}

coef_Killing_tp
coef_iniziali

iterazioni
lasso=lars(X1,Y1)
 plot(lasso)
 lasso$beta
 lasso[["lambda"]]
     
     
     
     
     # 
# 
# # PROVA DEL FOR CON K=1 
# 
# Xk_Prova <- X1[,3]                  # vettore del PRIMO regressore
# Xk_Prova
# Xk1_prova <- cbind(Xk_Prova)        # incolonno i valori di Xk => lo uso per il calcolo di Ck
# Xk1_prova
# XkRestanti_Prova <- X1[,-3]         # matrice dei restanti k-1 regressori
# XkRestanti_Prova
# betaRestanti_Prova <- coef_TOP[-3]  # vettore dei restanti k-1 parametri
# betaRestanti_Prova
# ak_Prova <- 2*(sum((Xk_Prova)^2))   # calcolo di ak
# ak_Prova
# 
# residui_Prova <- Y1 - (XkRestanti_Prova%*%betaRestanti_Prova) # faccio Y - X*beta (per i rimanenti k-1 regressori/parametri)
# residui_Prova
# Ck_Prova <- 2*(sum(residui_Prova * Xk1_prova))
# Ck_Prova
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
