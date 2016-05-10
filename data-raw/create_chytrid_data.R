chytrid <- read.csv("data-raw/coldinwarm.csv")
names(chytrid) <- c("time","count")
devtools::use_data(chytrid)
