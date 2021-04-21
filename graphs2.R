mu1 <- 50
sig1 <- 15

mu2 <- 20
sig2 <- 5

mu3 <- 40
sig3 <- 5

mu4 <- 50
sig4 <- 20

mu5 <- 10
sig5 <- 30

f1 <- function(x) {exp(-(x-mu1)^2/(2*sig1^2)) /(sig1*sqrt(2*pi)) }
f2 <- function(x) {exp(-(x-mu2)^2/(2*sig2^2)) /(sig2*sqrt(2*pi)) }
f3 <- function(x) {exp(-(x-mu3)^2/(2*sig3^2)) /(sig3*sqrt(2*pi)) }
f4 <- function(x) {exp(-(x-mu4)^2/(2*sig4^2)) /(sig4*sqrt(2*pi)) }
f5 <- function(x) {exp(-(x-mu5)^2/(2*sig5^2)) /(sig5*sqrt(2*pi)) }

f <- function(x) f1(x)
e <- function(x) 70*f5(x)

lwst <- .039
expan <- .00165
expit2 <- function(x) abs(1/(1+exp(-x))-.5+lwst)
se <- function(x) expit2(expan*(50-x)^2)
tau <- function(x) { 10*(f2(x-30)-f3(x-30)) }

r_e <- function(x) { 10*(f2(x)-f4(x)) }
i <- function(x) {.5*log((x+1)/(1-x)) }

p <- function(x) {gamma*(1-abs(r_e(x)))^i(alpha) }
#gam <- function(x){(1+exp(.1*(x-60)))^-1}

p1 <- function(x) {gam(x)*abs(tau(x))^i(alpha)*(1-se(x))^i(beta) }


asre_f <- function(x) { p(x)*f(x)*tau(x)*r_e(x) }
asre_f0 <- function(x) { p0(x)*f(x)*tau(x)*(ifelse(tau(x)>0,1,0)-e(x)) }
asre_f1 <- function(x) { p1(x)*f(x)*tau(x)*(ifelse(tau(x)>0,1,0)-e(x)) }

are_f1 <- function(x) { f(x)*tau(x)*(ifelse(tau(x)>0,1,0)-e(x)) }

p2 <- function(x) {gamma*abs(tau(x))^i(alpha) }
asre_f2 <- function(x) { p2(x)*f(x)*tau(x)*r_e(x) }

sim<- seq(10,100,by=.1)

wl <- 10
dev.new(width=wl, height=wl*3, pointsize=7)
par(mar=c(7,4,1,1))
par(mfrow=c(3,1))
par(xpd=FALSE)

# plot 1

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(-1,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))
abline(h=0)
points(sim, 40*sapply(sim, f), type="l", lwd=2, col="#374E55FF")
#points(sim, sapply(sim, r_e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, tau), type="l", lwd=2, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)+se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)-se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")

par(xpd=TRUE)
legend("topright", box.lty = 0, bg="white",
       legend=c("f(x) . 40", expression(hat(e)(x)), expression(hat(tau)(x))),
       lwd=2, col=c("#374E55FF", "#00A1D5FF","#DF8F44FF"))

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))

jamacol1 <- c("#79AF97FF", "#79C897")[1]
jamacol2 <- c("#B24745FF", "#C84745")[1]
jamacol <- c(jamacol2, jamacol1)

n <- 4
it <- 0
p0 <- function(x) {gamma }

for (k in seq(0,1,length=n)[-c(1,n)]) {
        gamma <- k
it <- it+1
lines(sim, sapply(sim, p0), lwd=2, col=jamacol[it])
}

legend(8,1.06, box.lty = 0, title=expression(p(x)==gamma),
       legend=rev(c(expression(paste(gamma, "=1/3 ")),
                expression(paste(gamma, "=2/3 ")))),
       lwd=2, col=rev(jamacol))


plot(1, type="n", bty="n", xaxt="n", yaxt="n", xlim=range(sim), ylim=c(0,.18), xlab="age", ylab="")
axis(side=1, at = seq(10,100, by=15))
axis(side=2, at = seq(0,.18, length=4), las=1)

n <- 4
it <- 0
asre_val <- list()

for (k in seq(0,1,length=n)[-c(1,n)]) {
        gamma <- k
        it <- it+1
        lines(sim, sapply(sim, function(x) integrate(asre_f0, lower=10, upper=x)$value), lwd=2, col=jamacol[it])
        asre_val[[it]] <- integrate(asre_f0, lower=10, upper=100)$value
}

par(xpd=FALSE)
asre_val <- sapply(asre_val,c)
are_val <- integrate(are_f1, lower=10, upper=100)$value

segments(0,are_val,100,are_val, lty=2)

legend(8,.13, box.lty = 0, title=expression(ASRE(r^stoch)), bg="white",
       legend=rev(format(round(asre_val, digits=2), nsmall = 2)), 
       lwd=2, col=rev(jamacol))
text(13.5, are_val+.01, paste("ARE(r) =", format(round(are_val, digits=2), nsmall = 2)) )

### plot 2

dev.new(width=wl, height=wl*3, pointsize=7)
par(mar=c(7,4,1,1))
par(mfrow=c(3,1))
par(xpd=FALSE)

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(-1,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))
abline(h=0)
points(sim, 40*sapply(sim, f), type="l", lwd=2, col="#374E55FF")
#points(sim, sapply(sim, r_e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, tau), type="l", lwd=2, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)+se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)-se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")


legend("topright", box.lty = 0, bg="white",
       legend=c("f(x) . 40", expression(hat(e)(x)), expression(hat(tau)(x))),
       lwd=2, col=c("#374E55FF", "#00A1D5FF","#DF8F44FF"))

par(xpd=TRUE)
plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))

jamacol1 <- c("#79AF97FF", "#79C897")[1]
jamacol2 <- c("#B24745FF", "#C84745")[1]
jamacol <- c(jamacol1, jamacol2)

it <- 0
p0 <- function(x) (1-se(x))^i(beta)

for (k in c(.80,.95)) {
                beta <- k
                it <- it+1
                lines(sim, sapply(sim, p0), lwd=2, col=jamacol[it])
        }

legend(8,1.03, box.lty = 0, title=expression(p(x)==(1-se[tau(x)])^legit(beta)),
       legend=c(expression(paste(beta, "=0.80")),
                expression(paste(beta, "=0.95"))),
       lwd=2, col=jamacol)


plot(1, type="n", bty="n", xaxt="n", yaxt="n", xlim=range(sim), ylim=c(0,.18), xlab="age", ylab="")
axis(side=1, at = seq(10,100, by=15))
axis(side=2, at = seq(0,.18, length=4), las=1)

n <- 5
it <- 0
asre_val <- list()

for (k in c(.80,.95)) {
                beta <- k
                it <- it+1
                lines(sim, sapply(sim, function(x) integrate(asre_f0, lower=10, upper=x)$value), lwd=2, col=jamacol[it])
                asre_val[[it]] <- integrate(asre_f0, lower=10, upper=100)$value
        }

asre_val <- sapply(asre_val,c)
are_val <- integrate(are_f1, lower=10, upper=100)$value

par(xpd=FALSE)
segments(0,are_val,100,are_val, lty=2)

legend(8,.13, box.lty = 0, title=expression(ASRE(r^stoch)), bg="white",
       legend=format(round(asre_val, digits=2), nsmall = 2), 
       lwd=2, col=jamacol)
text(13.5, are_val+.01, paste("ARE(r) =", format(round(are_val, digits=2), nsmall = 2)) )

### plot 3
dev.new(width=wl, height=wl*3, pointsize=7)
par(mar=c(7,4,1,1))
par(mfrow=c(3,1))
par(xpd=FALSE)

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(-1,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))
abline(h=0)
points(sim, 40*sapply(sim, f), type="l", lwd=2, col="#374E55FF")
#points(sim, sapply(sim, r_e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, tau), type="l", lwd=2, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)+se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)-se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")

par(xpd=TRUE)
legend("topright", box.lty = 0, bg="white",
       legend=c("f(x) . 40", expression(hat(e)(x)), expression(hat(tau)(x))),
       lwd=2, col=c("#374E55FF", "#00A1D5FF","#DF8F44FF"))

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))

jamacol1 <- c("#79AF97FF", "#79C897")
jamacol2 <- c("#B24745FF", "#C84745")
jamacol <- c(jamacol1, jamacol2)

it <- 0
p0 <- function(x) (1-se(x))^i(beta) * (-x/90+10/9)^i(alpha)

for (j in c(.80,.95)) {
        alpha <- j
for (k in c(.80,.95)) {
        beta <- k
        it <- it+1
        lines(sim, sapply(sim, p0), lwd=2, col=jamacol[it])
}
}

legend(70,1.03, box.lty = 0, title=expression(p(x)==(-x/90+10/9)^legit(alpha)~(1-se[tau(x)])^legit(beta)),
       legend=c(expression(paste(alpha, "=0.80 ", beta, "=0.80")),
                expression(paste(alpha, "=0.80 ", beta, "=0.95")),
                expression(paste(alpha, "=0.95 ", beta, "=0.80")),
                expression(paste(alpha, "=0.95 ", beta, "=0.95"))),
       lwd=2, col=jamacol)


plot(1, type="n", bty="n", xaxt="n", yaxt="n", xlim=range(sim), ylim=c(0,.18), xlab="age", ylab="")
axis(side=1, at = seq(10,100, by=15))
axis(side=2, at = seq(0,.18, length=4), las=1)

it <- 0
asre_val <- list()

for (j in c(.80,.95)) {
        alpha <- j
        for (k in c(.80,.95)) {
                beta <- k
        it <- it+1
        lines(sim, sapply(sim, function(x) integrate(asre_f0, lower=10, upper=x)$value), lwd=2, col=jamacol[it])
        asre_val[[it]] <- integrate(asre_f0, lower=10, upper=100)$value
        }
}

asre_val <- sapply(asre_val,c)
are_val <- integrate(are_f1, lower=10, upper=100)$value

par(xpd=FALSE)
segments(0,are_val,100,are_val, lty=2)

legend(8,.13, box.lty = 0, title=expression(ASRE(r^stoch)), bg="white",
       legend=format(round(asre_val, digits=2), nsmall = 2), 
       lwd=2, col=jamacol)
text(13.5, are_val+.01, paste("ARE(r) =", format(round(are_val, digits=2), nsmall = 2)) )

### plot 4
dev.new(width=wl, height=wl*3, pointsize=7)
par(mar=c(7,4,1,1))
par(mfrow=c(3,1))
par(xpd=FALSE)

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(-1,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))
abline(h=0)
points(sim, 40*sapply(sim, f), type="l", lwd=2, col="#374E55FF")
#points(sim, sapply(sim, r_e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, tau), type="l", lwd=2, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)+se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)-se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")

par(xpd=TRUE)
legend("topright", box.lty = 0, bg="white",
       legend=c("f(x) . 40", expression(hat(e)(x)), expression(hat(tau)(x))),
       lwd=2, col=c("#374E55FF", "#00A1D5FF","#DF8F44FF"))

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))

jamacol1 <- c("#79AF97FF", "#79C897")[1]
jamacol2 <- c("#B24745FF", "#C84745")[1]
jamacol <- c(jamacol1, jamacol2)

p0 <- function(x) as.numeric(sign((tau(x)+qt*se(x))*(tau(x)-qt*se(x)))==1)

it <- 0
for (j in c(0.45, .05) ) {
        qt <- qnorm(1-j/2)
                it <- it+1
                lines(sim, sapply(sim, p0)-(it-1)*.01, lwd=2, col=jamacol[it])
        }

legend(8,1.03, box.lty = 0, title=expression(p(x)==1[(tau(x)+q[1-alpha/2]*se(x))*(tau(x)-q[1-alpha/2]*se(x))>0]),
       legend=c(expression(paste(alpha, "=0.45 ")),
                expression(paste(alpha, "=0.05 "))),
       lwd=2, col=jamacol)


plot(1, type="n", bty="n", xaxt="n", yaxt="n", xlim=range(sim), ylim=c(0,.18), xlab="age", ylab="")
axis(side=1, at = seq(10,100, by=15))
axis(side=2, at = seq(0,.18, length=4), las=1)

it <- 0
asre_val <- list()

for (j in c(0.3, .05) ) {
        qt <- qnorm(1-j/2)
        it <- it+1
                lines(seq(10,100, by=.5), sapply(seq(10,100, by=.5)+.05, function(x) integrate(asre_f0, lower=10, upper=x)$value), lwd=2, col=jamacol[it])
                asre_val[[it]] <- integrate(asre_f0, lower=10, upper=100)$value
        }

asre_val <- sapply(asre_val,c)
are_val <- integrate(are_f1, lower=10, upper=100)$value

par(xpd=FALSE)
segments(0,are_val,100,are_val, lty=2)

legend(8,.13, box.lty = 0, title=expression(ASRE(r^stoch)), bg="white",
       legend=format(round(asre_val, digits=2), nsmall = 2), 
       lwd=2, col=jamacol)
text(13.5, are_val+.01, paste("ARE(r) =", format(round(are_val, digits=2), nsmall = 2)) )

### plot 5
dev.new(width=wl, height=wl*3, pointsize=7)
par(mar=c(7,4,1,1))
par(mfrow=c(3,1))
par(xpd=FALSE)

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(-1,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))
abline(h=0)
points(sim, 40*sapply(sim, f), type="l", lwd=2, col="#374E55FF")
#points(sim, sapply(sim, r_e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, e), type="l", lwd=2, col="#00A1D5FF")
points(sim, sapply(sim, tau), type="l", lwd=2, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)+se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")
points(sim, sapply(sim, function(x) tau(x)-se(x)), type="l", lwd=1, lty=3, col="#DF8F44FF")

par(xpd=TRUE)
legend("topright", box.lty = 0, bg="white",
       legend=c("f(x) . 40", expression(hat(e)(x)), expression(hat(tau)(x))),
       lwd=2, col=c("#374E55FF", "#00A1D5FF","#DF8F44FF"))

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))

jamacol1 <- c("#79AF97FF", "#79C897")[1]
jamacol2 <- c("#B24745FF", "#C84745")[1]
jamacol <- c(jamacol1, jamacol2)

p0 <- function(x) abs(1-ifelse(tau(x)>0,1,0)-e(x))^i(alpha)

it <- 0
for (j in c(1/3, 2/3) ) {
        alpha <- j
        it <- it+1
        lines(sim, sapply(sim, p0), lwd=2, col=jamacol[it])
}

legend(8,.6, box.lty = 0, title=expression(p(x)==abs(1-1[tau(x)>0]-e(x))^legit(alpha)),
       legend=c(expression(paste(alpha, "=1/3 ")),
                expression(paste(alpha, "=2/3 "))),
       lwd=2, col=jamacol)


plot(1, type="n", bty="n", xaxt="n", yaxt="n", xlim=range(sim), ylim=c(0,.18), xlab="age", ylab="")
axis(side=1, at = seq(10,100, by=15))
axis(side=2, at = seq(0,.18, length=4), las=1)

it <- 0
asre_val <- list()

for (j in c(1/3, 2/3) ) {
        alpha <- j
        it <- it+1
        lines(seq(10,100, by=.5), sapply(seq(10,100, by=.5), function(x) integrate(asre_f0, lower=10, upper=x)$value), lwd=2, col=jamacol[it])
        asre_val[[it]] <- integrate(asre_f0, lower=10, upper=100)$value
}

asre_val <- sapply(asre_val,c)
are_val <- integrate(are_f1, lower=10, upper=100)$value

par(xpd=FALSE)
segments(0,are_val,100,are_val, lty=2)

legend(8,.13, box.lty = 0, title=expression(ASRE(r^stoch)), bg="white",
       legend=format(round(asre_val, digits=2), nsmall = 2), 
       lwd=2, col=jamacol)
text(13.5, are_val+.01, paste("ARE(r) =", format(round(are_val, digits=2), nsmall = 2)) )
