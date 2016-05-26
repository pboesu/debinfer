#this function should combine the make.hypers and log.prior.params
setup.mcmc.params<-function(hyper=NULL, w.p=NULL, samp.type=NULL,
                            prop.var=NULL, prop.mean=NULL, start=NULL,
                            joint=FALSE, n.joint=NULL, scale=NULL){


  ## first load the parameters for the priors
  if(is.null(hyper)) hyper<-make.hypers()

  ## Need to change these bits if parameters are only being proposed
  ## individually. Additionally, most of these can be changed by
  ## passing in a new vector with the appropriate name and format as
  ## an arguement to the function. However, most bits for any sets of
  ## parameters that are being proposed jointly must be changed
  ## manually below.

  if(is.null(w.p)) w.p <- c("sd.Z", "sr", "fs", "ds", "eta",
                            "Tmin", "Tmax", "muz") ##alpha
  if(is.null(samp.type)) samp.type<-list(sd.Z="rw", alpha = "ind", sr="ind",
                                         fs="ind",
                                         ds="ind", eta="ind", Tmin="rw",
                                         Tmax="ind", muz="rw")
  if(is.null(prop.var)) prop.var<-list(sd.Z=c(4,5), alpha = 0.25, sr=0.0001,
                                       fs=0.001,
                                       ds=c(4,5), eta=2, Tmin=0.001,
                                       Tmax=1.5, muz=c(9,10))
  if(is.null(prop.mean)) prop.mean<-list(sd.Z=NULL, alpha = NULL, sr=NULL,
                                         fs=NULL,
                                         ds=NULL, eta=NULL, Tmin=NULL,
                                         Tmax=NULL, muz=NULL)
  if(is.null(start)) start<-c(sd.Z=0.18, alpha = 2, sr=0.2, fs=0.3,
                              ds=0.5, eta=19, Tmin=35,
                              Tmax=42, muz=0.24)


  ## The following lines must be changed if parameters are being
  ## proposed jointly

  ## for the oldline sims
  #  w1<-c("eta", "sr", "muz", "Tmin")
  #  cor1<-matrix(c(1, -0.2, 0.8, 0.4,
  #                 -0.2, 1, -0.5, 0.65,
  #                 0.8, -0.5, 1, 0.1,
  #                 0.4, 0.65, 0.1, 1), ncol=4, byrow=TRUE)
  #  cor1<-matrix(c( 1, -0.38,   0.76, 0.095,
  #                 -0.38, 1, -0.47, 0.69,
  #                 0.76, -0.47,  1, 0.14,
  #                 0.95, 0.69,  0.14, 1),  ncol=4, byrow=TRUE)
  #  hyp1<-c(hyper$eta, hyper$sr, hyper$muz, hyper$Tmin)
  #  m1<-c(3, 1, 0.5, 4)
  #  t1<-"rw"

  ##   for the newline sims
  w1<-c("sr", "Tmin")
  cor1<-matrix(c(1, 0.65,
                 0.65, 1), ncol=2, byrow=TRUE)
  hyp1<-c(hyper$sr, hyper$Tmin)
  m1<-c(0.8, 0.6)
  t1<-"rw"


  w2<-c("eta", "fs", "muz")
  cor2<-matrix(c(1, -0.57, 0.21,
                 -0.57, 1, 0.61,
                 0.21, 0.61, 1), ncol=3, byrow=TRUE)
  hyp2<-c(hyper$eta, hyper$fs, hyper$muz)
  m2<-c(0.8, 0.6, 0.5)
  t2<-"rw"

  w1<-c("sr", "Tmin","eta", "fs", "muz")
  cor1<-matrix(c(1, 0.763,  0.109, -0.0277, -0.154,
                 0.763, 1,  0.578, -0.0337,  0.406,
                 0.109, 0.578, 1, -0.393, 0.861,
                 -0.0277, -0.0337, -0.393,  1, 0.00163,
                 -0.154, 0.406,  0.861,  0.00163,  1), ncol=5, byrow=TRUE)
  hyp1<-c(hyper$sr, hyper$Tmin, hyper$eta, hyper$fs, hyper$muz)
  m1<-c(0.8, 0.6, 0.8, 0.6, 0.5)
  t1<-"rw"






  ws<-list(w1)#, w2)
  ms<-list(m1)#, m2)
  cors<-list(cor1)#, cor2)
  hyps<-list(hyp1)#, hyp2)
  types<-list(t1)#, t2)

  ## The code below here builds the appropriate structure with all the
  ## info needed by the mcmc to propose samples, etc.

  all<-list()

  if(joint){
    for(i in 1:n.joint){
      all[[i]] <- list(params = ws[[i]],
                       ##var = diag(x=c(0.00055, 0.0005), nrow=2),
                       mean   = ms[[i]],
                       var    = make.sigNN(var=prop.var[ws[[i]]],
                                           cor=cors[[i]], scale),
                       hyp    = matrix(hyps[[i]], nrow=length(ws[[i]]),
                                       byrow=TRUE),
                       start  = start[ws[[i]]],
                       type   = types[[i]]
      )
    }
  }

  if(!joint) n.joint<-0
  if(!is.null(w.p)){
    n.params <- n.joint+length(w.p)

    for(i in (n.joint+1):n.params){
      all[[i]]<-list(params=w.p[i-n.joint])
    }

    for(i in 1:n.params){
      pp<-all[[i]]$params
      if(length(pp)==1){
        if(is.element(pp, w.p)){
          all[[i]]$type<-samp.type[[pp]]
          all[[i]]$mean<-prop.mean[[pp]]
          all[[i]]$var<-prop.var[[pp]]
          all[[i]]$hyp<-hyper[[pp]]
          all[[i]]$start<-start[[pp]]
        }
        else{
          all[[i]]$var<-0
          ##all[[i]]$hyp<-0
          ##all[[i]]$
        }
      }
    }

  }
  return(all)
}




# > str(mcmc.p)
# List of 6
# $ :List of 5
# ..$ params: chr "fs"
# ..$ type  : chr "ind"
# ..$ var   : num 0.005
# ..$ hyp   : num [1:2] 1 1
# ..$ start : num 0.99
# $ :List of 5
# ..$ params: chr "ds"
# ..$ type  : chr "rw"
# ..$ var   : num [1:2] 3 4
# ..$ hyp   : num [1:2] 1 1
# ..$ start : num 1.71
# $ :List of 5
# ..$ params: chr "muz"
# ..$ type  : chr "rw"
# ..$ var   : num [1:2] 1 2
# ..$ hyp   : num [1:2] 5 1
# ..$ start : num 1.02
# $ :List of 5
# ..$ params: chr "eta"
# ..$ type  : chr "rw"
# ..$ var   : num [1:2] 4 5
# ..$ hyp   : num [1:2] 1 0.25
# ..$ start : num 36.4
# $ :List of 5
# ..$ params: chr "Tmin"
# ..$ type  : chr "rw"
# ..$ var   : num 0.1
# ..$ hyp   : num [1:2] 40 20
# ..$ start : num 3.9
# $ :List of 5
# ..$ params: chr "sr"
# ..$ type  : chr "ind"
# ..$ var   : num 1
# ..$ hyp   : num [1:2] 5 1
# ..$ start : num 5.11
