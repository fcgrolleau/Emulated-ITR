are_boot <- function(d, i=1:nrow(d)) {
        z<-d[i,]
        return(with(z, mean(Y*(r-E)*(A-E)/(E*(1-E))) ))
}

y1a0r0_boot <- function(d, i=1:nrow(d)) {
        z<-d[i,]
        return(with(z, 
                    sum(Y*A*r*(1-E)/E) / sum( (1-A)*r )) )
}

area0r0_boot <- function(d, i=1:nrow(d)) {
        z<-d[i,]
        return(with(z, sum(Y*A*r*(1-E)/E) / sum( (1-A)*r )) -
                       with(z[z$A==0 & z$r!=z$A,], mean(Y)) )
}

y1a1r0_boot <- function(d, i=1:nrow(d)) {
        z<-d[i,]
        return(with(z, 
                    sum(Y*(1-A)*E*(1-r)/(1-E)) / sum( (1-r)*A )
        ))
}

area1r0_boot <- function(d, i=1:nrow(d)) {
        z<-d[i,]
        return(with(z, 
                    sum(Y*(1-A)*E*(1-r)/(1-E)) / sum( (1-r)*A )) -
                       with(z[z$A==1 & z$r!=z$A,], mean(Y))
        )
}


