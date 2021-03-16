load("imputedpooleddata2.RData")
library(mice)
imp1 <- rbind(complete(imp_pooled_precoce,1), complete(imp_pooled_tardif,1))

## cross validation for random forest hyperparameter
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
     ylab="Cross validated AUC",
     xaxt="n")
axis(1, at = seq(1, nfeat, by = 1), las=1)

for (i in 1:nfeat){
        points(rep(i,nfolds), cv.disc[,i], col="grey", pch=19, cex=.3)
}
abline(v=which.max(apply(cv.disc, 2, mean)), col="green", lwd=3, lty=2)

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

# Compute the ARE estimator according to equation
are_hat <- with(imp1idealicu, 
     mean(Y*(r-E)*(A-E)/(E*(1-E)))
     )

## bootstrap
res <- boot(imp1idealicu, are_boot, R=999)
are_hat_ci <- boot.ci(res)$bca[4:5]

# print nice results
are_table <- paste0(round(c(are_hat, are_hat_ci)*100,2), "%")
names(are_table) <- c("Estimated ARE", "95% CI [ -", "- ]")
