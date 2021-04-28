#n simulated patients
nsim <- 2000

# set mean of variables simalar to those of the EXCEL trial
m <- c(66,75,57,0,0,0,0,0,0,0,20.5)
# name these variables
# NB:
# lmcd="Left main coronary artery disease only"
# tvd="Three-vessel disease only"
# if lmcd=F & tvd=F the patient both that is she has "left main coronary artery disease plus three-vessel disease"
names(m) <- c("age", "clearance", "lvef", "copd", "pvd", "dm", "insulin", "smoking", "tvd", "lmcd", "syntax")

# set variances of variables simalar to those of the EXCEL trial
var <- c(9.6^2, 17^2, 9.6^2,1,1,1,1,1,1,1,6.1^2)
names(var) <- c("age", "clearance", "lvef", "copd", "pvd", "dm", "insulin", "smoking", "tvd", "lmcd", "syntax")

# set proportion of categorical variables simalar to those of the EXCEL trial
prop <- c(7,9.5,29,7.7,22,60, 5)/100
names(prop) <- c("copd", "pvd", "dm", "insulin", "smoking", "tvd", "lmcd")

cov <- matrix(rep(0, length=length(m)^2), ncol = length(m))
colnames(cov) <- names(m)
rownames(cov) <- names(m)
diag(cov) <- rep(1, length(m))

cov["age","clearance"] <- -.7
cov["age","lvef"] <- -.5
cov["age","pvd"] <- .2
cov["age","tvd"] <- .3
cov["age","syntax"] <- .4

cov["smoking","clearance"] <- -.4
cov["smoking","lvef"] <- -.2
cov["smoking","copd"] <- .8
cov["smoking","pvd"] <- .7
cov["smoking","tvd"] <- .6
cov["smoking","syntax"] <- .7

cov["dm","clearance"] <- -.8
cov["dm","lvef"] <- -.4
cov["dm","insulin"] <- .8

sym <- function (x) {
        mat <- matrix(NA, ncol=ncol(x), nrow=nrow(x))
        colnames(mat) <- colnames(x)
        rownames(mat) <- rownames(x)
        for(j in 1:ncol(x)) {
        for(i in 1:nrow(x)) {
        if(x[i,j]==x[j,i]) {
                mat[i,j] <- x[i,j]
        } else {
                mat[i,j] <- ifelse(x[i,j]!=0, x[i,j], x[j,i])
        }
        }
        }
return(mat)
}

cov <- sym(cov)
#cov2 <- as.matrix(Matrix::nearPD(cov, corr=T, do2eigen=FALSE)$mat)

library(mvtnorm)

set.seed(6479)
simdat <- as.data.frame(rmvnorm(n=nsim, mean=rep(0, length(m)), sigma=cov))
colnames(simdat) <- names(m)

iter <- 0
simdat2 <- as.data.frame(sapply(simdat, function(x) {iter <<-iter+1
        x*sqrt(var[iter]) + m[iter] }
       ))

apply(simdat2, 2, mean)
apply(simdat2, 2, sd)

for (i in names(prop)){
simdat2[,i] <- simdat2[,i]<qnorm(prop[i])
}

# correct for insulin (no non-diabetic patients take insulin)
# proportion of diabetic patients on insulin in EXCEL
probs_insulin_given_dm <- 7.7/29
simdat2[,"insulin"] <- FALSE
simdat2[simdat2$dm==T,"insulin"] <- as.logical(rbinom(sum(simdat2$dm==T), size=1, prob = probs_insulin_given_dm))

apply(simdat2, 2, mean)
apply(simdat2, 2, sd)

apply(simdat2[,names(prop)],2, mean)*100
cor(simdat2)
summary(simdat2)

tableone::CreateTableOne(vars=names(m), data=simdat2)

expit <- function(x) 1/(1+exp(-x))
beta0 <- 1.3
betas <- c(-.6,-.3,1.1)
names(betas) <- c("age", "copd", "syntax")
simdat2$ps <- as.vector(with(simdat2, expit(beta0 + cbind(age/10, copd, syntax/10) %*% betas)))
summary(simdat2$ps)

set.seed(6547)
simdat2$cabg <- rbinom(nrow(simdat2), size=1, prob = simdat2$ps)

tableone::CreateTableOne(vars=names(m), strata="cabg", data=simdat2)


syntax2020 <- function(x, cabg) { x <- sapply(x,c)
        1-exp(-.243*exp(.99*(.72*x[["age"]]/10-0.07*ifelse(x[["clearance"]]>90,90,x[["clearance"]])/10-0.31*ifelse(x[["lvef"]]>50,50,x[["lvef"]])/10+.48*x[["copd"]]+.73*x[["pvd"]]+.2*x[["dm"]]+.46*x[["insulin"]]+.66*x[["smoking"]])-.4*cabg*x[["tvd"]]-.08*cabg*x[["lmcd"]]-.1*(1-cabg)*x[["lmcd"]]+.16*(1-cabg)*(x[["syntax"]]-29)/10-2.8))
}

# Reproduce case No 1 from figure 4 in the original paper
case1 <- data.frame(age=74, clearance = 38.6, lvef = 40, copd = 0, pvd = 0, dm = 0, insulin = 0, smoking = 1, tvd = 0, lmcd = 1, syntax = 11)
syntax2020(case1, cabg = 0)
syntax2020(case1, cabg = 1)

# Reproduce case No 2 from figure 4 in the original paper
case2 <- data.frame(age=59, clearance = 67.6, lvef = 67, copd = 0, pvd = 0, dm = 1, insulin = 0, smoking = 0, tvd = 1, lmcd = 0, syntax = 10)
syntax2020(case2, cabg = 0)
syntax2020(case2, cabg = 1)

# Reproduce case No 3 from figure 4 in the original paper
case3 <- data.frame(age=69, clearance = 72.5, lvef = 55, copd = 0, pvd = 0, dm = 1, insulin = 1, smoking = 0, tvd = 1, lmcd = 0, syntax = 50)
syntax2020(case3, cabg = 0)
syntax2020(case3, cabg = 1)

simdat2$preds_ttt <- apply(simdat2, 1, function(x) syntax2020(x, cabg = x["cabg"])) 
simdat2$preds_cabg <- apply(simdat2, 1, function(x) syntax2020(x, cabg = 1)) 
simdat2$preds_pci <- apply(simdat2, 1, function(x) syntax2020(x, cabg = 0))
simdat2$tau <- simdat2$preds_cabg-simdat2$preds_pci

tableone::CreateTableOne(strata="cabg", data=simdat2)

simdat2[sample(nsim, size = 1),]

### Compute and plot ASREs from this simulated dataset under different scenarios

# treatment (cabg) effect: function tau
tau <- function(x) {
        return(syntax2020(x, cabg = 1)-syntax2020(x, cabg = 0))
}

# propensity score : e(x) in the oracle case (ps and prognostic models known and unbiased)
e <- function(x) as.numeric(expit(beta0 + t(betas) %*% sapply(c(x["age"]/10, x["copd"], x["syntax"]/10),c) ))

# r(x)-e(x): function r_e
# NB by definition r(x) = I{tau_hat(x)>0}
# Here we assume tau_hat is unbiased and tau_hat(x)=tau(x)
r_e <- function(x) ifelse(tau(x)>0,1,0) - e(x)

# legit function i
legit <- function(x) {.5*log((x+1)/(1-x)) }

# cognitive biais: function p
p0 <- function(x) (1-abs(r_e(x)))^legit(alpha)

# implementation if the confidence interval doesn't cross 0, no implementation otherwise: function p1
p1 <- function(x) as.numeric(sign((tau(x)+qt*se(x))*(tau(x)-qt*se(x)))==1)

# summand in the asre integral: asre_score functions for each p function
asre_score0 <- function(x) p0(x)*tau(x)*r_e(x)
asre_score1 <- function(x) p1(x)*tau(x)*r_e(x)

# are_score function
are_score <- function(x) tau(x)*r_e(x)

# Compute ARE
are <- mean(apply(simdat2, 1, function(x) {x <- x[names(m)]
        are_score(x)}))

# Set precision for the lines in the coming plot
prec <- 100

# asre0 will save the ASRE for different values of cognitive biais alpha between 0 and 1
asre0 <- list()

# perc0 will save the percentage of patients implementing the ITR 
# for different values of cognitive biais alpha between 0 and 1
perc0 <- list()

# Compute asre0 and perc0 for different values of cognitive bias between 0 and 1
# with precision depending on prec
iter <- 0
for (j in seq(0,1,length=prec+1)) {
iter <- iter+1
alpha <- j
asre0[[iter]] <-mean(
                apply(simdat2, 1, function(x) {x <- x[names(m)]
        asre_score0(x)})
        )
perc0[[iter]] <-mean(
        apply(simdat2, 1, function(x) {x <- x[names(m)]
        p0(x)})
)
print(paste0(round(iter/(prec+1),2)*100,"%"))
}
asre0 <- sapply(asre0,c)

## Now let's simulate standard errors of the ITEs as the authors of the original paper do no provide a way to compute them
# Below I run PCA on patients characteristics
# and then give larger standard errors to "rare" patients (patient with high |pca1| )
# and to patients less likely to be randomized in a trial 
# (patients with PS far off from .5 as this means equipoise is less likely)

acp <- PCA(scale(simdat2[,c(1:11)]))
simdat2$s_pca1 <- as.numeric(scale(acp$ind$coord[,"Dim.1"]))

equipose_coef <- .2
outlier_coef <- .1
se <- function(x) { x <- x[c(names(m),"s_pca1")] 
                equipose_coef*(abs( e(x) - .5 )) + outlier_coef*abs(x["s_pca1"])
                }
simdat2$se <- apply(simdat2, 1, se)
summary(simdat2)

# examine the patients with the highest simulated standard error for the ITE
simdat2[which.max(simdat2$se),]
# this is a 31yo smoking patient with diabetes and peripheral vascular disease : quite uncommon in this population
# his PS is .83 very likely to have surgery rather than PCI and so less lekly to be randomized i a trial like the original syntax
# from with syntax 2020 is derived


# asre1 will save the ASRE for different values of alpha threshold between 0 and 1
asre1 <- list()

# perc1 will save the percentage of patients implementing the ITR 
# for different values of alpha threshold between 0 and 1
perc1 <- list()

# Again, prec is the precision for the lines in the coming plot
prec <- 100

# Let's now compute asre1 and perc1 for different values of alpha threshold between 0 and 1
# with precision depending on prec
iter <- 0
for (j in seq(0,1,length=prec+1) ) {
        qt <- qnorm(1-j/2)
        iter <- iter+1
        asre1[[iter]] <-mean(
                apply(simdat2, 1, asre_score1)
        )
        perc1[[iter]] <-mean(
                apply(simdat2, 1, p1)
        )
        print(paste0(round(iter/(prec+1),2)*100,"%"))
}
asre1 <- sapply(asre1, c)


### get quantitites for the plots
t <- 11
imp_txt0 <- paste0(round(sapply(perc0,c),2)*100,"%")[trunc(seq(1, prec+1, length=t))[-c(1, length(seq(1, prec+1, length=t)))]]
imp_txt1 <- paste0(round(sapply(perc1,c),2)*100,"%")[trunc(seq(1, prec+1, length=t))[-c(1, length(seq(1, prec+1, length=t)))]]
x_pos0 <- seq(0,1,length=prec+1)[trunc(seq(1, prec+1, length=t))[-c(1, length(seq(1, prec+1, length=t)))]]
x_pos1 <- x_pos0 
y_pos0 <- asre0[trunc(seq(1, prec+1, length=t))[-c(1, length(seq(1, prec+1, length=t)))]]
y_pos1 <- asre1[trunc(seq(1, prec+1, length=t))[-c(1, length(seq(1, prec+1, length=t)))]]

### plot the results
wl <- 5
dev.new(width=2*wl, height=wl, pointsize=7, noRStudioGD = TRUE)
par(mar=c(4,4,4,1))
par(mfrow=c(1,2))
par(xpd=FALSE)
jamacol <- c("#B24745FF", "#79AF97FF", "#374E55FF")

par("xpd"=T)
plot(asre0~seq(0,1,length=prec+1), type="l", col=jamacol[1], bty="n", las=1, xaxs = "i", yaxs = "i",
     xlim=c(0,1), ylim=c(0,are), yaxt="n",
     xlab=expression(alpha), ylab="ASRE", lwd=2)
axis(2, at=seq(0, are, length=5), labels=round(seq(0, are, length=5),3), las=2)
segments(0,are, 2, are,lty=3, col=jamacol[3])
points(x_pos0, y_pos0, pch=15, col="white", cex=3)
for (j in 1:(t+1)){
        text(x_pos0[j], y_pos0[j], imp_txt0[j], cex=1)
}
mtext(paste0("ARE = ", format(round(are,3), nsmall = 3)), side=3, line=.05, at=.1, cex=1)
title(expression(p(x)==(1-abs(1[tau(x)>0]-e(x)))^legit(alpha)))

plot(asre1~seq(0,1,length=prec+1), type="l", col=jamacol[2], bty="n", las=1, xaxs = "i", yaxs = "i",
     xlim=c(0,1), ylim=c(0,are), yaxt="n",
     xlab=expression(alpha), ylab="", lwd=2)
points(x_pos1, y_pos1, pch=15, col="white", cex=3)
segments(-3,are, 1, are,lty=3, col=jamacol[3])
for (j in 1:(t+1)){
        text(x_pos1[j], y_pos1[j], imp_txt1[j], cex=1)
}

title(expression(p(x)==1[(tau(x)+q[1-alpha/2]*se(x))~(tau(x)-q[1-alpha/2]*se(x))>0]))

par(mar=c(4,4,4,1))
par(mfrow=c(1,1))
dev.new(width=wl, height=wl, pointsize=7, noRStudioGD = TRUE)
plot(perc0, asre0, type="l", col=jamacol[1], bty="n", las=1, xaxs = "i", yaxs = "i",
     xlim=c(0,1), ylim=c(0,are), yaxt="n",
     xlab="No. of patients implementing the rule / No. of patients in the whole population", ylab="ASRE", lwd=2)
axis(2, at=seq(0, are, length=5), labels=round(seq(0, are, length=5),3), las=2)
lines(perc1, asre1, type="l", col=jamacol[2], lwd=2)
segments(0,are, 2, are,lty=3, col=jamacol[3])
mtext(paste0("ARE = ", format(round(are,3), nsmall = 3)), side=3, line=.05, at=.1, cex=1)

legend(.55, are/4, box.lty = 0, bg="white",
      legend=c(expression(p(x)==1[(tau(x)+q[1-alpha/2]*se(x))~(tau(x)-q[1-alpha/2]*se(x))>0]),
               expression(p(x)==(1-abs(1[tau(x)>0]-e(x)))^legit(alpha))),
      lwd=2, col=c(jamacol[2], jamacol[1]))
#dev.copy2pdf(file="fig7.0.pdf")
