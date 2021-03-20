load("imputedpooleddata2.RData")
library(mice)
imp1 <- rbind(complete(imp_pooled_precoce,1), complete(imp_pooled_tardif,1))

## k-fold cross validation to pick the random forest hyperparameter
library(caret)
library(pROC)
nfolds <- 10
nfeat <- 23
fold <- createFolds(imp1[imp1$etude=="akiki" & imp1$bras=="STRATEGIE PRECOCE",]$etat.censureJ60,
                    k=nfolds,
                    list = FALSE)

rf.akiki.early <- list()
tb <- list()
c.ext.er <- list()
cv.er <- list()
auc.er <- list()
cv.disc<- list()

for (m in 1:nfeat){  
for (kf in 1:nfolds){        
rf.akiki.early[[kf]] <- randomForest(I(factor(etat.censureJ60))~
                                            age+sexe+hypertension+diabetes+cirrhosis+
                                            respiratory_disease+cancer+aides+immunosupressive_drug+
                                            organ_graft+respi_sofa+hemodynamic_sofa+bilirubin_sofa+
                                            platelet_sofa+gcs_sofa+creat_b+sofa_e+weight_e+creat_e+
                                            urea_e+pot_e+bicar_e+ph_e,
                                    mtry=m, ntree=1000,
                                    importance =F,
                                    data=imp1[imp1$etude=="akiki" & imp1$bras=="STRATEGIE PRECOCE",][fold!=kf,])

y.early.ext <- factor(imp1[imp1$etude=="akiki" & imp1$bras=="STRATEGIE PRECOCE", "etat.censureJ60"][fold==kf])
yhat.early.ext <-  factor(predict(rf.akiki.early[[kf]], newdata=imp1[imp1$etude=="akiki" & imp1$bras=="STRATEGIE PRECOCE",][fold==kf,]))
yhat.early.ext.prob <- predict(rf.akiki.early[[kf]], newdata=imp1[imp1$etude=="akiki" & imp1$bras=="STRATEGIE PRECOCE",][fold==kf,],  type="prob")[,2]
tb[[kf]] <- table(y.early.ext, yhat.early.ext)
auc.er[[kf]] <- as.numeric(auc(y.early.ext, yhat.early.ext.prob))
c.ext.er[[kf]] <- sum(diag(tb[[kf]][nrow(tb[[kf]]):1, ]))/sum(tb[[kf]])
}
cv.er[[m]] <- unlist(c.ext.er)
cv.disc[[m]] <- unlist(auc.er)
print(paste0(round(m*kf/(nfolds*nfeat)*100,2), "%"))
}

par(mfrow=c(1,2))
cv.er <- sapply(cv.er, c)
plot(apply(cv.er, 2, mean), ylim = c(0,1), col= "red", pch=19,
     xlab="m",
     ylab="Misclassification Error",
     xaxt="n")
axis(1, at = seq(1, nfeat, by = 1), las=1)

for (i in 1:nfeat){
        points(rep(i,nfolds), cv.er[,i], col="grey", pch=19, cex=.3)
}
abline(v=which.min(apply(cv.er, 2, mean)), col="green", lwd=3, lty=2)

cv.disc <- sapply(cv.disc, c)
plot(apply(cv.disc, 2, mean), ylim = c(0,1), col= "red", pch=19,
     xlab="m",
     ylab="AUC",
     xaxt="n")
axis(1, at = seq(1, nfeat, by = 1), las=1)

for (i in 1:nfeat){
        points(rep(i,nfolds), cv.disc[,i], col="grey", pch=19, cex=.3)
}
abline(v=which.max(apply(cv.disc, 2, mean)), col="green", lwd=3, lty=2)

# OK it seems that lots if interaction between the variables doesn't help much
# here while it's likely to overfit the training data
# from now on we arbitrarily pick the hyperparameter as ~ sqrt(p) without crossvalidation

### Get single imputed data for AKIKI & IDEAL-ICU trials
library(dplyr)
imp1akiki <- imp1 %>% filter(etude=="akiki")
imp1idealicu <- imp1 %>% filter(etude=="idealicu")

### AKIKI X-learner
## Step 1
## Estimate the response functions
# here hyper parameter is ~ sqrt(p)=4
mu.akiki.early <- randomForest(I(factor(etat.censureJ60))~
                                             age+sexe+hypertension+diabetes+cirrhosis+
                                             respiratory_disease+cancer+aides+immunosupressive_drug+
                                             respi_sofa+hemodynamic_sofa+bilirubin_sofa+
                                             platelet_sofa+gcs_sofa+creat_b+weight_e+creat_e+
                                             urea_e+pot_e+bicar_e+ph_e,
                                     mtry=4, ntree=1000,
                                     importance =T,
                                     data=imp1akiki[imp1akiki$bras=="STRATEGIE PRECOCE",])

mu.akiki.late <- randomForest(I(factor(etat.censureJ60))~
                                       age+sexe+hypertension+diabetes+cirrhosis+
                                       respiratory_disease+cancer+aides+immunosupressive_drug+
                                       respi_sofa+hemodynamic_sofa+bilirubin_sofa+
                                       platelet_sofa+gcs_sofa+creat_b+weight_e+creat_e+
                                       urea_e+pot_e+bicar_e+ph_e,
                               mtry=4, ntree=1000,
                               importance =T,
                               data=imp1akiki[imp1akiki$bras=="STRATEGIE D ATTENTE",])
#diagnostics
importance(mu.akiki.early)
importance(mu.akiki.late)
varImpPlot(mu.akiki.early)
varImpPlot(mu.akiki.late)

## Step 2 
## Impute the treatment effects for the individuals in the treated group,
## based on the control-outcome estimator
imp1akiki$mupred <- NA
imp1akiki$mutepred <- NA
imp1akiki[imp1akiki$bras=="STRATEGIE PRECOCE",]$mupred <- predict(mu.akiki.late, newdata=imp1akiki[imp1akiki$bras=="STRATEGIE PRECOCE",],  type="prob")[,2]
imp1akiki[imp1akiki$bras=="STRATEGIE PRECOCE",]$mutepred <- with(imp1akiki[imp1akiki$bras=="STRATEGIE PRECOCE",], etat.censureJ60 - mupred)

## Impute the treatment effects for the individuals in the control group, 
## based on the treatment-outcome estimator
imp1akiki[imp1akiki$bras=="STRATEGIE D ATTENTE",]$mupred <- predict(mu.akiki.early, newdata=imp1akiki[imp1akiki$bras=="STRATEGIE D ATTENTE",],  type="prob")[,2]
imp1akiki[imp1akiki$bras=="STRATEGIE D ATTENTE",]$mutepred <- with(imp1akiki[imp1akiki$bras=="STRATEGIE D ATTENTE",], mupred - etat.censureJ60)

## Step 3
## Use Random Forests to estimate tau(x) in two ways
tau.akiki.early <- randomForest(mutepred~
                                       age+sexe+hypertension+diabetes+cirrhosis+
                                       respiratory_disease+cancer+aides+immunosupressive_drug+
                                       respi_sofa+hemodynamic_sofa+bilirubin_sofa+
                                       platelet_sofa+gcs_sofa+creat_b+weight_e+creat_e+
                                       urea_e+pot_e+bicar_e+ph_e,
                               mtry=4, ntree=1000,
                               importance =T,
                               data=imp1akiki[imp1akiki$bras=="STRATEGIE PRECOCE",])

tau.akiki.late <- randomForest(mutepred~
                                        age+sexe+hypertension+diabetes+cirrhosis+
                                        respiratory_disease+cancer+aides+immunosupressive_drug+
                                        respi_sofa+hemodynamic_sofa+bilirubin_sofa+
                                        platelet_sofa+gcs_sofa+creat_b+weight_e+creat_e+
                                        urea_e+pot_e+bicar_e+ph_e,
                                mtry=4, ntree=1000,
                                importance =T,
                                data=imp1akiki[imp1akiki$bras=="STRATEGIE D ATTENTE",])

#diagnostics
importance(tau.akiki.early)
importance(tau.akiki.late)
varImpPlot(tau.akiki.early)
varImpPlot(tau.akiki.late)

## Step 4
## Define the CATE estimate by a weighted average of the two estimates in step 2
imp1akiki$cate <- NA
w_ps_akiki <- mean(imp1akiki$bras=="STRATEGIE PRECOCE")

imp1akiki$cate <- w_ps_akiki*predict(tau.akiki.early, newdata = imp1akiki) +
        (1-w_ps_akiki)*predict(tau.akiki.late, newdata = imp1akiki)


### IDEAL-ICU X-learner
## Step 1
## Estimate the response functions
# here hyper parameter is ~ sqrt(p)=4
mu.idealicu.early <- randomForest(I(factor(etat.censureJ60))~
                                       age+sexe+hypertension+diabetes+cirrhosis+
                                       respiratory_disease+cancer+aides+immunosupressive_drug+
                                       respi_sofa+hemodynamic_sofa+bilirubin_sofa+
                                       platelet_sofa+gcs_sofa+creat_b+weight_e+creat_e+
                                       urea_e+pot_e+bicar_e+ph_e,
                               mtry=4, ntree=1000,
                               importance =T,
                               data=imp1idealicu[imp1idealicu$bras=="STRATEGIE PRECOCE",])

mu.idealicu.late <- randomForest(I(factor(etat.censureJ60))~
                                      age+sexe+hypertension+diabetes+cirrhosis+
                                      respiratory_disease+cancer+aides+immunosupressive_drug+
                                      respi_sofa+hemodynamic_sofa+bilirubin_sofa+
                                      platelet_sofa+gcs_sofa+creat_b+weight_e+creat_e+
                                      urea_e+pot_e+bicar_e+ph_e,
                              mtry=4, ntree=1000,
                              importance =T,
                              data=imp1idealicu[imp1idealicu$bras=="STRATEGIE D ATTENTE",])
#diagnostics
importance(mu.idealicu.early)
importance(mu.idealicu.late)
varImpPlot(mu.idealicu.early)
varImpPlot(mu.idealicu.late)

## Step 2 
## Impute the treatment effects for the individuals in the treated group,
## based on the control-outcome estimator
imp1idealicu$mupred <- NA
imp1idealicu$mutepred <- NA
imp1idealicu[imp1idealicu$bras=="STRATEGIE PRECOCE",]$mupred <- predict(mu.idealicu.late, newdata=imp1idealicu[imp1idealicu$bras=="STRATEGIE PRECOCE",],  type="prob")[,2]
imp1idealicu[imp1idealicu$bras=="STRATEGIE PRECOCE",]$mutepred <- with(imp1idealicu[imp1idealicu$bras=="STRATEGIE PRECOCE",], etat.censureJ60 - mupred)

## Impute the treatment effects for the individuals in the control group, 
## based on the treatment-outcome estimator
imp1idealicu[imp1idealicu$bras=="STRATEGIE D ATTENTE",]$mupred <- predict(mu.idealicu.early, newdata=imp1idealicu[imp1idealicu$bras=="STRATEGIE D ATTENTE",],  type="prob")[,2]
imp1idealicu[imp1idealicu$bras=="STRATEGIE D ATTENTE",]$mutepred <- with(imp1idealicu[imp1idealicu$bras=="STRATEGIE D ATTENTE",], mupred - etat.censureJ60)

## Step 3
## Use Random Forests to estimate tau(x) in two ways
tau.idealicu.early <- randomForest(mutepred~
                                        age+sexe+hypertension+diabetes+cirrhosis+
                                        respiratory_disease+cancer+aides+immunosupressive_drug+
                                        respi_sofa+hemodynamic_sofa+bilirubin_sofa+
                                        platelet_sofa+gcs_sofa+creat_b+weight_e+creat_e+
                                        urea_e+pot_e+bicar_e+ph_e,
                                mtry=4, ntree=1000,
                                importance =T,
                                data=imp1idealicu[imp1idealicu$bras=="STRATEGIE PRECOCE",])

tau.idealicu.late <- randomForest(mutepred~
                                       age+sexe+hypertension+diabetes+cirrhosis+
                                       respiratory_disease+cancer+aides+immunosupressive_drug+
                                       respi_sofa+hemodynamic_sofa+bilirubin_sofa+
                                       platelet_sofa+gcs_sofa+creat_b+weight_e+creat_e+
                                       urea_e+pot_e+bicar_e+ph_e,
                               mtry=4, ntree=1000,
                               importance =T,
                               data=imp1idealicu[imp1idealicu$bras=="STRATEGIE D ATTENTE",])

#diagnostics
importance(tau.idealicu.early)
importance(tau.idealicu.late)
varImpPlot(tau.idealicu.early)
varImpPlot(tau.idealicu.late)

## Step 4
## Define the CATE estimate by a weighted average of the two estimates in step 2
imp1idealicu$cate <- NA
w_ps_idealicu <- mean(imp1idealicu$bras=="STRATEGIE PRECOCE")

imp1idealicu$cate <- w_ps_idealicu*predict(tau.idealicu.early, newdata = imp1idealicu) +
        (1-w_ps_idealicu)*predict(tau.idealicu.late, newdata = imp1idealicu)


### Make IDEAL-ICU an observational study

## get prognosis predictions from both mu models
imp1idealicu$prognosis_pred <- apply(
cbind(predict(mu.idealicu.late, newdata=imp1idealicu,  type="prob")[,2],
predict(mu.idealicu.early, newdata=imp1idealicu,  type="prob")[,2]),
1,mean)

## Allocate strategy according to prognosis predictions with higher probability
## of allocation to early strategy when prognosis predictions is bad

imp1idealicu$reallocated_tt <- 
sapply(imp1idealicu$prognosis_pred, function(x) rbinom(n=1, size=1, prob=x))

imp1idealicu$reallocated_tt <- 
ifelse(imp1idealicu$reallocated_tt==1, "STRATEGIE PRECOCE", "STRATEGIE D ATTENTE")

## reallocate outcomes for the patients whose reallocated treatment differs from the original
imp1idealicu$reallocated_outcome <- NA
for (i in 1:nrow(imp1idealicu)) {
        if(with(imp1idealicu[i,], reallocated_tt==bras)) {
        imp1idealicu[i,]$reallocated_outcome <- imp1idealicu[i,]$etat.censureJ60
        } else {
                if(with(imp1idealicu[i,], reallocated_tt=="STRATEGIE PRECOCE")) {
                        imp1idealicu[i,]$reallocated_outcome <-
                                rbinom(n=1, size=1, prob=predict(mu.idealicu.early, newdata=imp1idealicu[i,],  type="prob")[,2])
                } else {
                        if(with(imp1idealicu[i,], reallocated_tt=="STRATEGIE D ATTENTE")) {
                                imp1idealicu[i,]$reallocated_outcome <-
                                        rbinom(n=1, size=1, prob=predict(mu.idealicu.late, newdata=imp1idealicu[i,],  type="prob")[,2])
                        
                }
        }
}
}

### Emulated trial estimation of ARE = E[Y(1)-Y(0)] in observational IDEAL-ICU
# Create relevant variables
imp1idealicu$A <- ifelse(imp1idealicu$reallocated_tt=="STRATEGIE PRECOCE", 1, 0)
imp1idealicu$Y <-imp1idealicu$reallocated_outcome

# Using the SOFA as only variable in a slightly mispecified propensity score model
ps_model <- glm(Y~sofa_e, data=imp1idealicu, family = "binomial")
imp1idealicu$E <-predict(ps_model, type = "response")

# rule from AKIKI X-learner : little r as follows
# if CATE is negative recommend early strategy (treatment 1) else late (treatment O)
imp1idealicu$cate.akiki <- w_ps_akiki*predict(tau.akiki.early, newdata = imp1idealicu) +
        (1-w_ps_akiki)*predict(tau.akiki.late, newdata = imp1idealicu)
imp1idealicu$r <-  ifelse(imp1idealicu$cate.akiki<0,1,0)

# add R as defined in the manuscript
imp1idealicu$R <- ifelse(with(imp1idealicu, A==r), 1, 0)

# Compute the ARE estimator according to equation
are_hat <- with(imp1idealicu, 
     mean(Y*(r-E)*(A-E)/(E*(1-E)))
     )

## bootstrap
res <- boot(imp1idealicu, are_boot, R=999)
are_hat_ci <- boot.ci(res)$bca[4:5]

# print results in a table
are_table <- paste0(round(c(are_hat, are_hat_ci)*100,2), "%")
names(are_table) <- c("Estimated ARE", "95% CI [ -", "- ]")
are_table

### Emulated trial estimation of E[Y(1) | A = 0, R = 0] in observational IDEAL-ICU
y1a0r0_hat <- with(imp1idealicu, 
                   sum(Y*A*r*(1-E)/E) / sum( (1-A)*r )
)

## bootstrap
res <- boot(imp1idealicu, y1a0r0_boot, R=999)
y1a0r0_ci <- boot.ci(res)$bca[4:5]

### Emulated trial estimation of E[Y(1) - Y | A = 0, R = 0] in observational IDEAL-ICU
area0r0_hat <- with(imp1idealicu, sum(Y*A*r*(1-E)/E) / sum( (1-A)*r )) -
        with(imp1idealicu[imp1idealicu$A==0 & imp1idealicu$r!=imp1idealicu$A,], mean(Y))

res <- boot(imp1idealicu, area0r0_boot, R=999)
area0r0_ci <- boot.ci(res)$bca[4:5]
        
### Emulated trial estimation of E[Y(1) | A = 1, R = 0] in observational IDEAL-ICU
y1a1r0_hat <- with(imp1idealicu, 
                   sum(Y*(1-A)*E*(1-r)/(1-E)) / sum( (1-r)*A )
)

## bootstrap
res <- boot(imp1idealicu, y1a1r0_boot, R=999)
y1a1r0_ci <- boot.ci(res)$bca[4:5]

### Emulated trial estimation of E[Y(1) - Y | A = 1, R = 0] in observational IDEAL-ICU
area1r0_hat <- with(imp1idealicu, 
                   sum(Y*(1-A)*E*(1-r)/(1-E)) / sum( (1-r)*A )) -
        with(imp1idealicu[imp1idealicu$A==1 & imp1idealicu$r!=imp1idealicu$A,], mean(Y))

## bootstrap
res <- boot(imp1idealicu, area1r0_boot, R=999)
area1r0_ci <- boot.ci(res)$bca[4:5]


## Check if estimations are consistent with the equation 

with(imp1idealicu, mean(R==0))* (
        with(imp1idealicu[imp1idealicu$R==0,], mean(A==0))*area0r0_hat
                +
        with(imp1idealicu[imp1idealicu$R==0,], mean(A==1))*area1r0_hat
)

are_hat        

# Both are equal so there should be not error in the implementation
# Here the benefit of the ITR is mostly (/totally) driven by area1r0_hat

### Compute the ARE estimator from AIPW estimator (Tsiatis et al. eq 3.30 page 62)
# get propensity for receiving treatment according to the rule (eq 3.21 page 59)
imp1idealicu$Pi_d <- with(imp1idealicu, (E^r)*((1-E)^(1-r)) )

# Note Cd,i in the book corresponds to Ri under our notation 

# Compute are_aipw_hat with Pi_d from slightly mispesified propensity score (SOFA)
# and prognosis_pred as prognostic score model (we don't know if it's well specified 
# as we don't know the real prognostic score model here)

are_aipw_hat <- with(imp1idealicu,
     mean(
          R*Y/Pi_d -
          prognosis_pred*(R-Pi_d)/Pi_d        
     ) - 
     mean(Y)
     )

# Note that when the prognosis score model is 0 are_aipw_hat == are_hat
are_aipw_hat <- with(imp1idealicu,
                    mean(
                        R*Y/Pi_d -
                        0*(R-Pi_d)/Pi_d        
                     ) - 
                    mean(Y)
)

dec <- vector()
for (i in 1:20){
dec[i] <- round(are_aipw_hat,i)==round(are_hat,i)
}
length(which(dec))

# aipw_hat == are_hat to 16 decimal places here


# boot strap the original are_aipw_hat
are_aipw_hat <- with(imp1idealicu,
                     mean(
                        R*Y/Pi_d -
                        prognosis_pred*(R-Pi_d)/Pi_d        
                     ) - 
                     mean(Y)
)

res <- boot(imp1idealicu, are_aipw_boot, R=999)
are_aipw_ci <- boot.ci(res)$bca[4:5]


### Simulations for effect of stochastic implementation of ITR on Delta aipw

set.seed(56489)
nimp <- 100
nsim <- 100

are_s_aipw_sim <- matrix(NA, nrow=nimp, ncol=nsim)

row <- 0
for (i in seq(0,1,length=nimp)) {
stoch_p <- i
row <- row + 1
for (j in 1:nsim) {
        
imp1idealicu$P <- rbinom(nrow(imp1idealicu),1,stoch_p)
# create new stochastic rule
imp1idealicu$r_s <- apply(imp1idealicu[,c('r','A','P')], 1, function(x) x['r']^rbinom(1,1,x['P'])*x['A']^(1-rbinom(1,1,x['P'])) )

imp1idealicu$R_s <- ifelse(with(imp1idealicu, A==r_s), 1, 0)

# this equation is critical for the simulation see notes/manuscript
imp1idealicu$Pi_d_s <- with(imp1idealicu, ((E^r_s)*((1-E)^(1-r_s)))^P )

are_s_aipw_sim[row,j] <- 
with(imp1idealicu,
                     mean(
                             R_s*Y/Pi_d_s -
                             prognosis_pred*(R_s-Pi_d_s)/Pi_d_s        
                     ) - 
                     mean(Y) )

}
print(paste0(100*row*j/(nimp*nsim),"%"))
}
are_s_aipw_sim

## bootstrap

ci_mat <- matrix(NA, nrow=nimp, ncol=2)
row <- 0
for (i in seq(0,1,length=nimp)) {
        row <- row + 1
        temp <- boot(imp1idealicu, are_s_aipw_boot, R=49, nsim=2, stoch_p=i)
        ci_mat[row,] <- temp$t0+c(-1,1)*qnorm(.975)*sd(temp$t)
        print(paste0(100*row/nimp,"%"))
}

## plot the stochastic implementation of the rule with confidence intervals
dev.new(width=10, height=10, unit="in")
plot(seq(0,1,length=nimp), apply(are_s_aipw_sim, 1, mean), ylim=c(-.08,.02),
     main="Benefit For ITR Implementation in The IDEAL-ICU Population",
     xlab="Stochastic Rule Implementation",
     ylab=expression(paste(Delta, "aipw For Mortality at Day 60")), type="n", bty="n", xaxt='n', las=2)
axis(1, seq(0,1,by=.2), paste0(seq(0,1,by=.2)*100, "%"))
lines(seq(0,1,length=nimp), apply(are_s_aipw_sim, 1, mean), pch=19, lwd=2, col=rgb(0, 161, 213, maxColorValue=255))

polygon(c(seq(0,1,length=nimp), rev(seq(0,1,length=nimp))), c(ci_mat[,1], rev(ci_mat[,2])),
        col= rgb(0, 161, 213, maxColorValue=255, alpha=255*.6), border=NA)
abline(h=0, lwd=1, lty=2)
#dev.copy2pdf(file="stochastic_plot.pdf")


### Simulations for effect of stochastic implementation of ITR on DeltaA0R0 / DeltaA1R0

set.seed(9543)
nimp <- 10
nsim <- 1000

delta_ratio_sim <- matrix(NA, nrow=nimp, ncol=nsim)

row <- 0
for (i in seq(0,1,length=nimp)) {
        stoch_p <- i
        row <- row + 1
for (j in 1:nsim) {

imp1idealicu$P <- rbinom(nrow(imp1idealicu),1,stoch_p)
# create new stochastic rule
imp1idealicu$r_s <- apply(imp1idealicu[,c('r','A','P')], 1, function(x) x['r']^rbinom(1,1,x['P'])*x['A']^(1-rbinom(1,1,x['P'])) )
# create new R variable for the stochastic rule                
imp1idealicu$R_s <- ifelse(with(imp1idealicu, A==r_s), 1, 0)

area0r0_temp <- with(imp1idealicu, sum(Y*A*r_s*(1-E)/E) / sum( (1-A)*r_s )) -
        with(imp1idealicu[imp1idealicu$A==0 & imp1idealicu$r!=imp1idealicu$A,], mean(Y))

area1r0_temp <- with(imp1idealicu, 
                    sum(Y*(1-A)*E*(1-r_s)/(1-E)) / sum( (1-r_s)*A )) -
        with(imp1idealicu[imp1idealicu$A==1 & imp1idealicu$r_s!=imp1idealicu$A,], mean(Y))
temp <- area0r0_temp/area1r0_temp
delta_ratio_sim[row, j] <- temp
}
print(paste0(100*j*row/(nimp*nsim), "%"))
}


###

# estimation of area0r0_temp and area1r0_temp are impossible in all the simulations 
# where the stochastic rule dictates that all treatment should be the same as in 
# standard of care
apply(delta_ratio_sim, 1, function(x) sum(is.nan(x)))

# we therefore ignore these simulations
apply(delta_ratio_sim, 1, mean, na.rm=T)

## plot the stochastic implementation of the rule with confidence intervals
dev.new(width=10, height=10, unit="in")
plot(seq(0,1,length=nimp)[-1], apply(delta_ratio_sim, 1, mean, na.rm=T)[-1],
     xlim=c(0,1), ylim=c(-.15, 1.15),
     main="What Drives The Benefit of ITR Implementation \nin The IDEAL-ICU Population",
     xlab="Stochastic Rule Implementation",
     ylab=expression(paste(Delta, "A0R0 / ", Delta, "A1R0")), type="n", bty="n", xaxt='n', las=2)
axis(1, seq(0,1,by=.2), paste0(seq(0,1,by=.2)*100, "%"))
points(seq(0,1,length=nimp)[-1], apply(delta_ratio_sim, 1, mean, na.rm=T)[-1], pch=19, col=rgb(0, 161, 213, maxColorValue=255))
lines(seq(0,1,length=nimp)[-1], apply(delta_ratio_sim, 1, mean, na.rm=T)[-1], pch=19, lwd=1, col=rgb(0, 161, 213, maxColorValue=255))
abline(h=1, lwd=1, lty=1)
abline(h=0, lwd=1, lty=2)
text(.3, .8, "For 100% implementation", pos=4)
text(.3, .75, expression(paste(Delta, "A0R0 =  0.004")), pos=4)
text(.3, .70, expression(paste(Delta, "A1R0 = -0.287")), pos=4)
text(.3, .50, "Benefit is only driven from implementing \nthe rule in the treated patients!", pos=4)
#dev.copy2pdf(file="stochastic_delta_ratio_plot.pdf")

