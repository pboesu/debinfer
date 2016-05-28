##traceplot for multiple chains
#obsolete when using plot.mcmc.list
plot_chains <- function(chains, nrow, ncol, cols = c('orange','red','darkgreen','cornflowerblue')){
  nch = length(chains)
  np = ncol(chains[[1]]$samps)
  if(nrow*ncol < np) stop("more parameters than plot can fit")
  w.p = names(chains[[1]]$samps)
  ranges = apply(do.call('rbind', (do.call('rbind', chains))), 2, range)
  par(mfrow = c(nrow, ncol))
  for (p in 1:np){
    plot(chains[[1]]$samps[,p], col = cols[1], main = w.p[p], type = 'l', ylim = ranges[,p])
    for (i in 2:nch){
      lines(chains[[i]]$samps[,p], col = cols[i], main = w.p[i])
    }
  }
}
