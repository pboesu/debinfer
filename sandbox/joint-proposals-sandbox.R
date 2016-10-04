#sandbox for joint proposals

debinfer_cov <- function(var.names, sigma=diag(length(names)), joint.block ){
  if(!is.character(var.names)) stop("var.names must be a character vector")
  if(class(sigma)!= "matrix" || !is.numeric(sigma)) stop("sigma must be a numeric matrix")
  if(any(dim(sigma)!=length(var.names))) stop("length(var.names) does not match dimensions of sigma")
  colnames(sigma)<-var.names
  rownames(sigma)<-var.names
  structure(list(sigma=sigma, joint.block = joint.block), class="debinfer_cov")
}

sigma1 <- debinfer_cov(c("r", "K"), sigma = matrix(c(1,0.25,0.25,1),2), "sigma1")


dmvnorm(x=c(0,0))
dmvnorm(x=c(0,0), mean=c(1,1))

sigma <- matrix(c(4,2,2,3), ncol=2)
x <- mvtnorm::rmvnorm(n=500, mean=c(0,0), sigma=sigma1$sigma)
xcolMeans(x)
var(x)

x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma, method="chol")
colMeans(x)
var(x)

plot(x)
