

### pairwise plot of marginal densities
### code following a stackexchange example

#' Pairwise posterior marginals
#'
#' Plots pairwise correlations of posterior marginals
#'
#' @param x a deBInfer_result object
#' @param trend logical, add loess smooth
#' @param scatter logical, add scatterplot of posterior samples
#' @param burnin integer, number of samples to discard from start of chain before plotting
#' @param medians logical, plot marginal medians on contour plot
#' @param ... further arguments to plot.default (the call that draws the scatter/contour plot)
#' @import MASS
#' @import RColorBrewer
#' @importFrom graphics abline contour hist layout lines par plot plot.default points text
#' @import stats
#' @export
pairs.debinfer_result <- function(x, trend = FALSE, scatter = FALSE, burnin=NULL, medians=TRUE, ...){
  if(!is.null(burnin)) x$samples <- window(x$samples, burnin, nrow(x$samples))
  np = ncol(x$samples)
  cors<-round(cor(x$samples),2) #correlations

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
    hist(x$samples[,i],main=colnames(x$samples)[i])

  # Write correlations to upper diagonal part of the graph
  # Again, tweak accordingly
  for(i in 1:(np-1))
    for(j in (i+1):np){
      plot.default(-1:1,-1:1, type = "n",xlab="",ylab="",xaxt="n",yaxt="n")
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
        plot.default(x$samples[,i],x$samples[,j], pch=16, cex=0.3, col='darkgrey', ...)
      } else {
        plot.default(range(x$samples[,i]), range(x$samples[,j]), type = 'n', ...)
      }
      z <- MASS::kde2d(x$samples[,i],x$samples[,j], n=20)
      contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
      if (medians) abline(h=median(x$samples[,j]), v=median(x$samples[,i]), lwd=2, lty = 2)
      if (trend == TRUE) {
        rows <- sample(nrow(x$samples), ceiling(0.05*nrow(x$samples)))
        loess_fit <- loess(x$samples[rows,j] ~ x$samples[rows,i])
        points(x$samples[rows,i], predict(loess_fit), col = "blue", pch=16, cex=0.5)
      }

    }
  # restore old par
  par(old.par)
  }

#' Plot posterior marginals and corresponding priors
#'
#' Plots posterior densities and the densities of the corresponding priors. The prior density is automatically evaluated for the range given by the x-axis limits of the plot (which defaults to the posterior support).
#'
#' @param result a deBInfer_result object
#' @param burnin numeric, number of samples to discard before plotting
#' @param param character, name of parameter to plot. "all" (default) plots all parameters
#' @param prior.col character color for prior density
#' @param n, integer, number of points at which to evaluate the prior density.
#' @param ... further arguments to coda::densplot
#' @import coda
#' @importFrom graphics abline contour hist layout lines par plot plot.default points text
#' @importFrom grDevices n2mfrow
#' @export
post_prior_densplot <- function(result, param="all", burnin=NULL, prior.col="red", n=1000, ...){

  #remove burnin if supplied
  if(!is.null(burnin)) result$samples <- window(result$samples, burnin, nrow(result$samples))
  if (param=="all"){
    # store old par
    old.par <- par(no.readonly=TRUE)
    #get number of parameters
    nplots <- ncol(result$samples)
    par(mfrow=n2mfrow(nplots))
    for (i in colnames(result$samples)){
      #plot posterior density
      coda::densplot(result$samples[,i],..., main = i)

      #construct prior densities
      dprior <- paste("d", result$all.params[[i]]$prior, sep="")
      #get x range form plot
      plot.range <- seq(par("usr")[1], par("usr")[2], length.out = n)
      #evaluate prior
      prior.dens <- do.call(dprior, c(list(x=plot.range), result$all.params[[i]]$hypers))
      #plot
      lines(plot.range, prior.dens, col=prior.col)
    }
    # restore old par
    par(old.par)
  } else {
    if (any(colnames(result$samples) == param )){
      i = param
      #plot posterior density
      coda::densplot(result$samples[,i],...)
      #construct prior densities
      dprior <- paste("d", result$all.params[[i]]$prior, sep="")
      #get x range form plot
      plot.range <- seq(par("usr")[1], par("usr")[2], length.out = n)
      #evaluate prior
      prior.dens <- do.call(dprior, c(list(x=plot.range), result$all.params[[i]]$hypers))
      #plot
      lines(plot.range, prior.dens, col=prior.col)
    } else {
      stop(paste(param, "is not a valid parameter name for this posterior sample."))
    }
  }
}

#' Plot inference outputs
#'
#' Plots the inference results from a debinfer_result object
#'
#' @param x a deBInfer_result object
#' @param plot.type character, which type of plot. Options are "coda" for coda::plot.mcmc, "post_prior" for deBInfer::post_prior_densplot.
#' @param ... further arguments to methods
#' @seealso \code{\link{post_prior_densplot}}, \code{\link[coda]{plot.mcmc}}, \code{\link{pairs.debinfer_result}}
#' @import coda
#' @export
plot.debinfer_result <- function(x, plot.type="coda", ...){
  if (plot.type=="coda") plot(x$samples, ...)
  if (plot.type=="post_prior") post_prior_densplot(x, ...)
  }

#' Plot posterior trajectory
#'
#' Plots the inference results from a debinfer_result object
#'
#' @param x a post_sim or post_sim_list object
#' @param plot.type character, which type of plot. Options are "ensemble" and "medianHDI".
#' @param col color, for plot.type="medianHDI" the first element is used for the median, the second for the HDI
#' @param ... further arguments to methods
#' 
#' @export
plot.post_sim_list <- function(x, plot.type="medianHDI", col = c("red","darkgrey"), lty = c(1,2), ...){
  if (plot.type=="medianHDI") {
    for (i in length(x$median)){
      plot(x$time, x$median[[i]], type='l', col = col[1], lty = lty[1], ylim = c(min(x$HDI[[i]][,1]),max(x$HDI[[i]][,2])), ylab = i, xlab = "time", ...)
      lines(x$time, x$HDI[[i]][,1], col=col[2], lty = lty[2])
      lines(x$time, x$HDI[[i]][,2], col=col[2], lty = lty[2])
    }
    }
  if (plot.type=="ensemble"){
    if(!inherits(x$sims[[1]], "deSolve")) stop("method not implemented for non deSolve solvers")
    do.call("plot", c(list(x = x$sims[[1]]),x$sims[2:length(x$sims)],list(col=col[1], lty=lty[1],...)))
  }
}
