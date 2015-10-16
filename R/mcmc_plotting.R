##traceplot for multiple chains

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

pretty_pairs <- function(samples, trend = FALSE, scatter = FALSE){
  np = ncol(samples)
  cors<-round(cor(samples),2) #correlations

  # make layout for plot layout
  laymat<-diag(1:np) #histograms
  laymat[lower.tri(laymat)]<-(np + 1):(np + (np^2 - np) / 2) #correlations
  laymat[upper.tri(laymat)]<-(np + (np^2 - np) / 2 + 1):np^2 #heatmaps

  layout(laymat) #define layout using laymat

  par(mar=c(2,2,2,2)) #define marginals etc.

  # Draw histograms, tweak arguments of hist to make nicer figures
  for(i in 1:np)
    hist(samples[,i],main=names(samples)[i])

  # Write correlations to upper diagonal part of the graph
  # Again, tweak accordingly
  for(i in 1:(np-1))
    for(j in (i+1):np){
      plot(-1:1,-1:1, type = "n",xlab="",ylab="",xaxt="n",yaxt="n")
      text(x=0,y=0,labels=paste(cors[i,j]),cex=2)
    }

  # Plot heatmaps, here I use kde2d function for density estimation
  # image function for generating heatmaps
  library(MASS)
  ## some pretty colors
  library(RColorBrewer)
  k <- 7
  my.cols <- (brewer.pal(9, "Blues")[3:9])

  ## compute 2D kernel density, see MASS book, pp. 130-131


  for(i in 2:np)
    for(j in 1:(i-1)){
      if (scatter == TRUE){
        plot(samples[,i],samples[,j], pch=16, cex=0.3, col='darkgrey')
      } else {
        plot(range(samples[,i]), range(samples[,j]), type = 'n')
      }
      z <- kde2d(samples[,i],samples[,j], n=20)
      contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
      abline(h=mean(samples[,j]), v=mean(samples[,i]), lwd=2, lty = 2)
      if (trend == TRUE) {
        rows <- sample(nrow(samples), ceiling(0.05*nrow(samples)))
        loess_fit <- loess(samples[rows,j] ~ samples[rows,i])
        points(samples[rows,i], predict(loess_fit), col = "blue", pch=16, cex=0.5)
      }

    }
  }


