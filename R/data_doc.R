#' Chytrid fungus data set
#'
#' Replicated spore counts of an experimental culture of the chytrid fungus \emph{Batrachochytrium dendrobatidis}.
#' This dataset is a subset of the observations from the experimental study conducted by Voyles et al. (2012).
#'
#' @name chytrid
#' @docType data
#' @format A data.frame with 76 rows and two columns
#' \describe{
#'   \item{time}{days since the start of the experiment}
#'   \item{count}{count of zoospores (x 1e4)}
#' }
#' @references Voyles et al. 2012, Ecol Evol 9:2241-2249 \doi{10.1002\%2Fece3.334}
#' @keywords data
NULL

#' Logistic growth data set
#'
#' Simulated data from the logistic growth model with N_0=0.1, r=0.1 and K=10
#'
#' @name logistic
#' @docType data
#' @format A data.frame with 36 rows and 3 columns
#' \describe{
#'   \item{time}{time since start of the model}
#'   \item{N_true}{Numerical solution of N_t }
#'   \item{N_noisy}{N_t with the addition of log-normal noise, where sdlog = 0.05}
#' }
#' @keywords data
NULL
