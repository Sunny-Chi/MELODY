#' Simulate Logistic Data with Covariates
#'
#' Generates logistic regression data incorporating covariates' effects, allowing for complex relationships between predictors and the outcome variable.
#'
#' @param sample_size Number of observations to generate.
#' @param mediator_size Number of mediators in the model.
#' @param p Total number of predictors.
#' @param a Coefficients affecting the relationship between X and M.
#' @param r Coefficient for the primary predictor X.
#' @param b Coefficients for mediators.
#' @param alpha1 Intercept in the logistic regression equation.
#' @param g Coefficient for additional predictors.
#' @param l Coefficient modeling the effect of Z on M.
#' @return A list containing the simulated dataset, mediator matrix, and probabilities.
#' @examples
#' sim_data_cov <- sim_logistic_cov(sample_size = 100, mediator_size = 2, p = 5)
#' @export

sim_logistic_cov= function(sample_size = sample_size, mediator_size = mediator_size, p = p, a = a, r = r, b = b, alpha1 = alpha1,g=g,l=l){
  X = rnorm(n = sample_size)           # generate x
  Z <-rnorm(n = sample_size)
  M<-matrix(NA,sample_size,mediator_size)
  for(i in 1:mediator_size){
    M[,i]<- a*X+l*Z+rnorm(sample_size)
  }

  beta = c(r,g, b)
  beta <- matrix(beta, ncol = 1)
  XX <- scale(cbind(X,Z,M))
  pr = 1/(1+exp(-(alpha1 + XX %*% beta)))  # generate y by equation 2
  y = rbinom(n = sample_size, size = 1, prob = pr)
  # generate noise mediators
  noise = matrix(rnorm(sample_size*(p - mediator_size)),nrow=sample_size,ncol=p - mediator_size)
  MT=M
  M=cbind(M,noise)
  colnames(M)=paste0('V', 1:ncol(M))
  dat <- data.frame(y, X, Z,M)
  return(list(dat=dat,MT=MT, pr=pr))
}


