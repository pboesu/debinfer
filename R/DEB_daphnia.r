

## Daphnia model (see Nisbet et al 2010 PTRS). Starvation rule should be completed

DEB.daphnia<-function(t,y,p){

   J.EAm <- p["J.EAm"]
   y.EX <- p["y.EX"]
   v <- p["v"]
   k.M <- p["k.M"]
   g <- p["g"]
   kap <- p["kap"]
   k.J <- p["k.J"]
   M.HP <- p["M.HP"]
   shape <- p["shape"]
   gamma <- p["gamma"]
   X.h <- p["X.h"]
   vol <- p["vol"]

   phi<-PHI(X.h, y[1]) #y[1]/(y[1]+X.h)

   F<-food(t,X.h,y[1])

   ## here Delta = -Delta in PTRS paper (ie, if >0 sufficient reserve
   ## to pay somatic maintenance)
   Delta <- - kap * g * J.EAm / (v * y[3]) * (( k.M * g * y[3]-y[2]*v)/ (y[2] + g))
   ##if(t%%1==0) print(paste("Delta = ", Delta))

   ## p.c, or nearly
   p.c<-(1 - kap) * y[3]^2 * y[2] * g * J.EAm / v * ( (v + k.M* y[3]) / (y[2] + g))

   if(Delta > 0){
     dy1 <- -J.EAm / (vol * y.EX) * y[3] * phi + F #dX
     dy2 <- v / y[3] * (phi-y[2]) #de
     dy3 <- 1/3 * ((y[2] * v - k.M * g * y[3])/(y[2] + g))#dL

     ## Development and reproduction
     if(y[4] < M.HP){
       dy4 <-  p.c  - k.J * y[4] #dMH
       dy5 <- 0 # dMR
     }
     else{
       dy4 <- 0
       dy5 <-  p.c - k.J * y[4]
     }
   }
   else{
     ##if(Delta>0){## J.EC > J.EJ + J.ES
       dy1 <- -J.EAm / (vol * y.EX) * y[3] * phi + F #dX
       dy2 <- v / y[3] * (phi-y[2]) #de
       dy3 <- 0 #dL
       ## Development and reproduction
       if(y[4] < M.HP){
         dy4 <- p.c - k.J * y[4] #dMH
             dy5 <- 0 # dMR
       }
       else{
         dy4 <- 0
         dy5 <- p.c - k.J * y[4]
       }
     ##}
     ##else{# death
     ##  dy1 <- F
     ##  dy2<- dy3 <- dy4 <- dy5 <-0
     ##}
   }


   ##print(paste("animal is dead at time ",t,sep=""))

   list(c(dy1,dy2,dy3,dy4,dy5))
 }


PHI<-function(X.h, x){
  p<-x/(x+X.h)
  return(p)
}


food<-function(t, X.h, x){

  if(t>1) f<-(1-x/(x+X.h))*(2-cos(t/10))##/2
  ##if(t%%7==0) f<-1000
  else f<-0

  ##f<-10-2*x
  return(f)
}



## Here are bits to set the parameters and initial conditions, and
## solve the system of odes, etc. There are also some plotting
## functions


setparams.DEB<-function(J.EAm = 0.00361, y.EX = 0.5, v = 1.48,
                        k.M = 0.105, g = 0.498, kap = 0.1,
                        k.J = 0.1, M.HP = 0.00575,
                        shape = 0.721, gamma = 630, X.h = 0.16,
                        vol = 0.1){

    params <- c(J.EAm = J.EAm, y.EX = y.EX, v = v, k.M = k.M, g = g,
                kap = kap, k.J = k.J*k.M, M.HP = M.HP, shape = shape,
                gamma = gamma, X.h = X.h, vol = vol)

  return(params)
}




#' Title
#' for this version of the model the initial enery density should be
#' set to the "parental" density, i.e. to food(1, X.h, X0), but for
#' the MCMC the "default is for e0==NULL, to indicate that e0 is
#' effectively also "proposed by the algorithm.
#'
#' @param X0
#' @param e0
#' @param L0
#' @param M.H0
#' @param M.R0
#'
#' @return initial states
#' @export
#'
#' @examples examples
setinits.DEB<-function(X0=0.2, e0=-Inf, L0 = 0.69, M.H0 = 0, M.R0 = 0){

  inits<-c(X = X0, e = e0, L = L0, M.H = M.H0, M.R = M.R0)
  return(inits)
}




## This function takes data from the forward simulator (either full,
## or some subset), and applies the observation model and noise to
## generate "data" that would be like the observed data. The
## "observations" are then: L=alpha*y3 + eps.L (eps.L ~N(0, sds$L),
## and Negg=gamma*y5 + eps.L (eps.L ~N(0, sds$Negg), where y3 is
## length and y5 is reproduction.

## The name of this function should remain the same, since other
## functions elsewhere call it.

add.noise<-function(data, sds, params){
  alpha<-params["shape"]
  gamma<-params["gamma"]
  ##print(c(alpha, gamma))
  t<-data[,'time']
  ##print(t)
  L<-rnorm(t, data[,'L']*alpha, sd=sds$L)
  Negg<-rnorm(t, data[,'M.R']*gamma, sd=sds$Negg)

  w<-which(L<0)
  L[w]<-0
  w<-which(Negg<0)
  Negg[w]<-0

  return(list(t=t, L=L, Negg=Negg))

}



