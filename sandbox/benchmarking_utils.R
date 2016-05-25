#benchmarking utilities
#devtools::install_github("GShotwell/convertr")
#library(convertr)
#explore_units()


plot_benchmark <- function(input, time.units="s"){
  levels(input$expr) <- ifelse(grepl("de_mcmc_rev", levels(input$expr)), "de_mcmc_rev", "de_mcmc")
  boxplot(convertr::convert(time, "ns", time.units)~expr, data= input, ylab = paste("Time (", time.units, ")", sep=""))
}
