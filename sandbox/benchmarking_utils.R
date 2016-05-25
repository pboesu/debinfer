#benchmarking utilities
#devtools::install_github("GShotwell/convertr")
#library(convertr)
#explore_units()


plot_benchmark <- function(input, time.units="s", ...){
  #rename levels
  levels(input$expr) <- ifelse(grepl("de_mcmc_rev", levels(input$expr)), "de_mcmc_rev", "de_mcmc")
  #convert time units
  input$time <- convertr::convert(input$time, "ns", time.units)
  #calculate plotting range
  #if (!exists("ylim")){
  #  if(yax.zero) {
  #    ylim = c(0, max(input$time))
  #  } else {
  #    ylim = c(min(input$time), max(input$time))
  #  }
  #}
  boxplot(time ~expr, data= input, ylab = paste("Time (", time.units, ")", sep=""), ...)
  grid()
}
