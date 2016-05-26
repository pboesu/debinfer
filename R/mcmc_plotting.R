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

### pairwise plot of marginal densities
### code following a stackexchange example
## check out the function in rethinking

#' pretty_pairs
#'
#' Plots pairwise correlations of posterior marginals
#'
#' @param result a deBInfer_result object
#' @param trend logical, add loess smooth
#' @param scatter logical, add scatterplot of posterior samples
#' @import MASS
#' @import RColorBrewer
#' @export
pairs.debinfer_result <- function(result, trend = FALSE, scatter = FALSE, burnin=NULL){
  if(!is.null(burnin)) result$samples <- window(result$samples, burnin, nrow(result$samples))
  np = ncol(result$samples)
  cors<-round(cor(result$samples),2) #correlations

  # store old par
  old.par <- par(no.readonly=TRUE)

  # make layout for plot layout
  laymat<-diag(1:np) #histograms
  laymat[lower.tri(laymat)]<-(np + 1):(np + (np^2 - np) / 2) #correlations
  laymat[upper.tri(laymat)]<-(np + (np^2 - np) / 2 + 1):np^2 #heatmaps

  layout(laymat) #define layout using laymat

  par(mar=c(2,2,2,2)) #define marginals etc.

  # Draw histograms, tweak arguments of hist to make nicer figures
  for(i in 1:np)
    hist(result$samples[,i],main=names(result$samples)[i])

  # Write correlations to upper diagonal part of the graph
  # Again, tweak accordingly
  for(i in 1:(np-1))
    for(j in (i+1):np){
      plot(-1:1,-1:1, type = "n",xlab="",ylab="",xaxt="n",yaxt="n")
      text(x=0,y=0,labels=paste(cors[i,j]),cex=2)
    }

  # Plot heatmaps, here I use kde2d function for density estimation
  # image function for generating heatmaps
  #library(MASS)
  ## some pretty colors
  #library(RColorBrewer)
  k <- 7
  my.cols <- (RColorBrewer::brewer.pal(9, "Blues")[3:9])

  ## compute 2D kernel density, see MASS book, pp. 130-131


  for(i in 2:np)
    for(j in 1:(i-1)){
      if (scatter == TRUE){
        plot.default(result$samples[,i],result$samples[,j], pch=16, cex=0.3, col='darkgrey')
      } else {
        plot.default(range(result$samples[,i]), range(result$samples[,j]), type = 'n')
      }
      z <- MASS::kde2d(result$samples[,i],result$samples[,j], n=20)
      contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
      abline(h=mean(result$samples[,j]), v=mean(result$samples[,i]), lwd=2, lty = 2)
      if (trend == TRUE) {
        rows <- sample(nrow(result$samples), ceiling(0.05*nrow(result$samples)))
        loess_fit <- loess(result$samples[rows,j] ~ result$samples[rows,i])
        points(result$samples[rows,i], predict(loess_fit), col = "blue", pch=16, cex=0.5)
      }

    }
  # restore old par
  par(old.par)
  }


