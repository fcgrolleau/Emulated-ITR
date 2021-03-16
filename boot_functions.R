are_boot <- function(d, i=1:nrow(d)) {
        z<-d[i,]
        return(with(z, mean(Y*(r-E)*(A-E)/(E*(1-E))) ))
}
