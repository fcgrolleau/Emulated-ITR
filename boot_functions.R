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

are_aipw_boot <- function(d, i=1:nrow(d)) {
        z<-d[i,]
        return(with(z,
                    mean(
                        R*Y/Pi_d -
                        prognosis_pred*(R-Pi_d)/Pi_d        
                    ) - 
                    mean(Y)
        ))
}


are_s_aipw_boot <- function(d, i=1:nrow(d), nsim=100, stoch_p=1) {
        z<-d[i,]
        temp <- vector()
        for (j in 1:nsim) {
                
                z$P <- rbinom(nrow(z),1,stoch_p)
                # create new stochastic rule
                z$r_s <- with(z, r^rbinom(1,1,P)*A^(1-rbinom(1,1,P)) )
                table(with(z, r==r_s))
                
                z$R_s <- ifelse(with(z, A==r_s), 1, 0)
                
                # this equation is critical for the simulation see notes/manuscript
                z$Pi_d_s <- with(z, ((E^r_s)*((1-E)^(1-r_s)))^P )
                
                temp[j] <- 
                        with(z,
                             mean(
                                R_s*Y/Pi_d_s -
                                prognosis_pred*(R_s-Pi_d_s)/Pi_d_s        
                             ) - 
                            mean(Y) ) 
                
        }
        return(mean(temp))
}
