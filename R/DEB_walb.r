

## WALB model (see Teixeira et al. 2014 J Se Res). Juvenile phase only

DEB.walb<-function(t,y,p){
  #"rename" variables. really there should be consistent naming for state and parameter vectors in the final version
  states <- y
  params <- p
  with(as.list(c(states, params)), {

    E_m = p_Am/ v;     #% J/cm^3, reserve capacity [E_m]
    g = E_G/ kap/ E_m; #% -, energy investment ratio

    TC = exp(T_A/T_ref - T_A/T_b);
    vT = v * TC;
    kT_J = k_J * TC;
    pT_Am = p_Am * TC;

    e = vT * E/ L^3/ pT_Am;             #% -, scaled reserve density;
    rT = vT * (e/ L - 1/ L_m)/ (e + g); #% 1/d, spec growth rate
    pT_C = E * (vT/ L - rT);            #% J/d, scaled mobilisation

    #[f_n] = get_f(t, f_slope, f_intercept);


    #derivatives
    dH = (1 - kap) * pT_C - kT_J * H; #% J
    dE = pT_Am * f_n * L^2 - pT_C;    #% J
    dL = rT * L/3;                    #% cm
    df_n = f_slope               #% -  #different from matlab implementation

    #return derivatives
    list(c(dH, dE, dL, df_n))
  })
}





## Here are bits to set the parameters and initial conditions, and
## solve the system of odes, etc. There are also some plotting
## functions


setparams.DEB.walb<-function(L_m = 26.62883424,
                        p_Am = 725.266318774218,
                        v = 0.029349082,
                        k_J = 2.7486143e-07,
                        kap = 0.9991795,
                        T_A = 15000,
                        T_ref = 293,
                        T_b = 312.6,
                        E_G = 13072.691,
                        f_slope = -0.0019625396,
                        f_intercept = 1.0659367){

    params <- c(L_m = L_m,
                p_Am = p_Am,
                v = v,
                k_J = k_J,
                kap = kap,
                T_A = T_A,
                T_ref = T_ref,
                T_b = T_b,
                E_G = E_G,
                f_slope = f_slope,
                f_intercept = f_intercept
                )

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
#'
#'
#' @examples examples
setinits.DEB.walb<-function(E_h = 2.756519348e+06,
                            L_h = 5.185484695e+00,
                            E_Hh = 2.360057500e+03,
                            f_n = 1.065936700e+00){

  inits<-c(H = E_Hh, E = E_h, L = L_h, f_n = f_n)
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

add.noise<-function(data, sds, params, ode.pars){
  alpha<-params["shape"]
  gamma<-params["gamma"]
  ##print(c(alpha, gamma))
  t<-data[,'time']
  ##print(t)
  L<-rnorm(t, data[,'L']*alpha, sd=sds$L)

  #wet weight. this is all hard coded now, should not be!
  w_E = 23.9 # molecular weight of reserve g mol^-1
  d_v = 0.5 # specific density of structure
  mu_E = 550000 # chemical potential of reserve J / mol
  wdratio = -1.37e-3 * data[, 'time'] + 2.09
  omega = unname(ode.pars['p_Am'] * w_E / (ode.pars['v'] * d_v * mu_E))
  Ww = data[,'L']^3 * (1 + data[,'f_n'] * omega) * d_v * wdratio
  w<-which(L<0)
  L[w]<-0
  w<-which(Ww<0)
  Ww[w]<-0

  return(list(t=t, L=L, Ww=Ww))

}



