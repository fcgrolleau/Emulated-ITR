mu1 <- 50
sig1 <- 15

mu2 <- 20
sig2 <- 5

mu3 <- 40
sig3 <- 5

mu4 <- 50
sig4 <- 20

f1 <- function(x) {exp(-(x-mu1)^2/(2*sig1^2)) /(sig1*sqrt(2*pi)) }
f2 <- function(x) {exp(-(x-mu2)^2/(2*sig2^2)) /(sig2*sqrt(2*pi)) }
f3 <- function(x) {exp(-(x-mu3)^2/(2*sig3^2)) /(sig3*sqrt(2*pi)) }
f4 <- function(x) {exp(-(x-mu4)^2/(2*sig4^2)) /(sig4*sqrt(2*pi)) }

f <- function(x) f1(x)       
tau <- function(x) { 10*(f2(x)-f3(x)) }
#tau <- function(x) { 10*(f2(x-30)-f3(x-30)) }
r_e <- function(x) { 10*(f2(x)-f4(x)) }
i <- function(x) {.5*log((x+1)/(1-x)) }
p <- function(x) {gamma*(1-abs(r_e(x)))^i(alpha) }
asre_f <- function(x) { p(x)*f(x)*tau(x)*r_e(x) }

sim<- seq(10,100,by=.1)

wl <- 10
dev.new(width=wl, height=wl*3, pointsize=7)
par(mar=c(7,4,1,1))
par(mfrow=c(3,1))
plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(-1,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))
abline(h=0)
points(sim, 40*sapply(sim, f), type="l", lwd=2, col="#374E55FF")
points(sim, sapply(sim, tau), type="l", lwd=2, col="#DF8F44FF")
points(sim, sapply(sim, r_e), type="l", lwd=2, col="#00A1D5FF")
legend("topright", box.lty = 0,
       legend=c("f(x) . 40", expression(tau(x)), "r(x)-e(x)"),
       lwd=2, col=c("#374E55FF", "#DF8F44FF", "#00A1D5FF"))

plot(1, type="n", bty="n", xaxt="n", xlim=range(sim), ylim=c(0,1), xlab="", ylab="", las=1)
axis(side=1, at = seq(10,100, by=15))

jamacol1 <- c("#79AF97FF", "#79C897")
jamacol2 <- c("#B24745FF", "#C84745")

n <- 4
gamma <- 1
it <- 0
for (j in seq(0,1,length=n)[-c(1,n)]) {
alpha <- j
it <- it+1
lines(sim, sapply(sim, p), lwd=2, col=jamacol1[it])
}

gamma <- .5
it <- 0
for (j in seq(0,1,length=n)[-c(1,n)]) {
        alpha <- j
        it <- it+1
        lines(sim, sapply(sim, p), lwd=2, col=jamacol2[it])
}

legend("topright", box.lty = 0, title="p(x)", bg="white",
       legend=c(expression(paste(gamma, "=1   ", alpha, "=1/3")),
                expression(paste(gamma, "=1   ", alpha, "=2/3")),
                expression(paste(gamma, "=1/2 ", alpha, "=1/3")),
                expression(paste(gamma, "=1/2 ", alpha, "=2/3"))),
       lwd=2, col=c(jamacol1, jamacol2))

plot(1, type="n", bty="n", xaxt="n", yaxt="n", xlim=range(sim), ylim=c(0,.003), xlab="age", ylab="")
#plot(1, type="n", bty="n", xaxt="n", yaxt="n", xlim=range(sim), ylim=c(-.004,.002), xlab="age", ylab="")
#abline(h=0)
axis(side=1, at = seq(10,100, by=15))
axis(side=2, at = seq(0,.003, by=.001), las=1)
#axis(side=2, at = seq(-.004,.002, by=.001), las=1)

gamma <- 1
n <- 4
it <- 0
for (j in seq(0,1,length=n)[-c(1,n)]) {
        alpha <- j
        it <- it+1
        lines(sim, sapply(sim, asre_f), lwd=2, col=jamacol1[it])
}

gamma <- .5
it <- 0
for (j in seq(0,1,length=n)[-c(1,n)]) {
        alpha <- j
        it <- it+1
        lines(sim, sapply(sim, asre_f), lwd=2, col=jamacol2[it])
}

asre <- list()
it <- 0
for (k in c(1,.5)) {
        gamma <- k
        for (j in seq(0,1,length=n)[-c(1,n)]) {
                alpha <- j
                it <- it+1
                asre[[it]] <- integrate(asre_f, lower=10, upper=100)$value
        }
}
asre <- unlist(asre)

legend("topright", box.lty = 0, title="ASRE", bg="white",
       legend=c(bquote(.(round(asre[1],2))~gamma==1~alpha==1/3),
                bquote(.(round(asre[2],2))~gamma==1~alpha==2/3),
                bquote(.(round(asre[3],2))~gamma==1/2~alpha==1/3),
                bquote(.(round(asre[4],2))~gamma==1/2~alpha==2/3),
                expression(paste(" ")) ),
       lwd=2, col=c(jamacol1, jamacol2, "white"))
#dev.copy2pdf(file="plottest_2.pdf")

asre <- list()
n <- 50
it <- 0
for (k in seq(0,1,length=n)) {
        gamma <- k
        for (j in seq(0,1,length=n)) {
                alpha <- j
                it <- it+1
                temp <- c(k,j, integrate(asre_f, lower=10, upper=100)$value)
                names(temp) <- c("gamma", "alpha", "asre")
                asre[[it]] <- temp
        }
}
asre<- as.data.frame(t(sapply(asre,c)))

library(akima)
library(plotly)
s <-  interp(x = asre$alpha, y = asre$gamma, z = asre$asre)
plot <- plot_ly(x = s$x, y = s$y, z = s$z) %>% add_surface()
plot
