#dde_solver_sandbox

voylespar <- c(sr = 1.2, fs = 0.95, ds =1.4, muz = 0.7, eta = 18, Tmin = 2.45)
voylesinit <- c(C=120, S=0, Z=0)
voylestimes <- 0:100

CSZ.dde<-function(t,y,p){

  sr    <-p["sr"]
  fs    <-p["fs"]
  ds    <-p["ds"]
  eta   <-p["eta"]
  Tmin  <-p["Tmin"]
  ##Tmax  <-p["Tmax"]
  muz <-p["muz"]

  Rs<-Ms<-0
  lag1<-lag2<-0

  if (t>Tmin){
    lag1<-pastvalue(t-Tmin)
    Rs <- sr*fs*lag1[1]
  }
  ##if (t>(Tmin+Tmax)){
  ##  lag2<-pastvalue(t-(Tmin+Tmax))
  ##  Ms <- sr*fs*exp(-ds*Tmax)*lag2[1]
  ##}

  phiZ <- eta*y[2]
  dy1 <- -(muz+sr)*y[1]
  dy2 <- Rs - Ms - ds*y[2]
  dy3 <- phiZ - (muz+sr)*y[3]

  if(y[1]<0) dy1<-0
  if(y[2]<0){
    dy2 <- Rs - Ms
    dy3 <- -(muz+sr)*y[3]
  }
  if(y[3]<0){
    dy3 <- dy3+(muz+sr)*y[3]
  }

  list(c(dy1,dy2,dy3))
}

CSZ.dede<-function(t,y,p){

  sr    <-p["sr"]
  fs    <-p["fs"]
  ds    <-p["ds"]
  eta   <-p["eta"]
  Tmin  <-p["Tmin"]
  ##Tmax  <-p["Tmax"]
  muz <-p["muz"]

  Rs<-Ms<-0
  lag1<-lag2<-0

  if (t>Tmin){
    lag1<-lagvalue(t-Tmin)
    Rs <- sr*fs*lag1[1]
  }
  ##if (t>(Tmin+Tmax)){
  ##  lag2<-pastvalue(t-(Tmin+Tmax))
  ##  Ms <- sr*fs*exp(-ds*Tmax)*lag2[1]
  ##}

  phiZ <- eta*y[2]
  dy1 <- -(muz+sr)*y[1]
  dy2 <- Rs - Ms - ds*y[2]
  dy3 <- phiZ - (muz+sr)*y[3]

  if(y[1]<0) dy1<-0
  if(y[2]<0){
    dy2 <- Rs - Ms
    dy3 <- -(muz+sr)*y[3]
  }
  if(y[3]<0){
    dy3 <- dy3+(muz+sr)*y[3]
  }

  list(c(dy1,dy2,dy3))
}


ddesoln <- PBSddesolve::dde(y = voylesinit, times = voylestimes, func = CSZ.dde, parms = voylespar)
deSolve:::plot.deSolve(ddesoln)

dedesoln <- dede(y = voylesinit, times = voylestimes, func = CSZ.dede, parms = voylespar)
plot(dedesoln)

for (i in 2:4){
  plot(ddesoln[,i], dedesoln[,i])
  abline(0,1, col='red')
}

#compare solver speed
microbenchmark(
  ddesoln <- dde(y = voylesinit, times = voylestimes, func = CSZ.dde, parms = voylespar),
  dedesoln <- as.data.frame(dede(y = voylesinit, times = voylestimes, func = CSZ.dede, parms = voylespar)),
  times=20
)

