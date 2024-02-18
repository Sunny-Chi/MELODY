
#' Simulate Logistic Data
#'
#' This function simulates logistic regression data. It supports simulation with or without covariates and allows adjustment of various simulation parameters.
#'
#' @param sample_size Integer, the number of samples to generate.
#' @param mediator_size Integer, the number of mediator variables.
#' @param p Integer, the total number of predictors.
#' @param a Numeric vector, coefficients for mediators.
#' @param r Numeric, coefficient for the primary predictor.
#' @param b Numeric vector, coefficients for the mediators.
#' @param alpha1 Numeric, intercept.
#' @param cov Optional covariance matrix.
#' @param g Additional parameters for `sim_logistic_cov`.
#' @param l Additional parameters for `sim_logistic_cov`.
#' @param nm1 Optional, additional parameters for specifying extra noise mediators.
#' @return A list containing the simulated dataset and additional simulation parameters.
#' @examples
#' sim_data <- sim_logistic_data(sample_size = 100, mediator_size = 2, p = 10)
#' @export
#' @import MASS DescTools SIS
sim_logistic_data = function(sample_size = 2000, mediator_size = 2, p = 1000, a = rep(0,2), r = 0, b = rep(0,2), alpha1 = 0,cov=NULL,g=g,l=l,nm1=NULL) {

  if(is.null(cov)==FALSE){
    return(sim_logistic_cov(sample_size = sample_size, mediator_size = mediator_size, p = p, a = a, r = r, b = b, alpha1 = alpha1,g=g,l=l))
  }else {
    if(is.null(nm1)){
      if(mediator_size==1){
        X <- rnorm(sample_size)
        noise = matrix(rnorm(sample_size*(p - mediator_size)),nrow=sample_size,ncol=p - mediator_size)
        MT <- rnorm(sample_size, mean=a*X, sd=a)
        M=cbind(MT,noise)
        colnames(M)=paste0('V', 1:ncol(M))
        beta = c(r, b)
        beta <- matrix(beta, ncol = 1)
        pr = 1/(1+exp(-(alpha1 + r*X + b*MT)))  # generate y by equation 2
        y = rbinom(n = sample_size, size = 1, prob = pr)
        dat <- data.frame(y, X, M)
        return(list(dat=dat,MT=MT, pr=pr))
      }else if(mediator_size==0){
        X <- rnorm(sample_size)
        noise = matrix(rnorm(sample_size*p),nrow=sample_size,ncol=p)

        M=noise
        colnames(M)=paste0('V', 1:ncol(M))
        pr = 1/(1+exp(-(alpha1 + r*X)))  # generate y by equation 2
        y = rbinom(n = sample_size, size = 1, prob = pr)
        dat <- data.frame(y, X, M)
        return(list(dat=dat,MT=NULL, pr=pr))
      }else{
        x = rnorm(n = sample_size)           # generate x
        #specify covariance matrix R among x and m
        R <- matrix(0, ncol = mediator_size,nrow = mediator_size)
        for(i in 1:mediator_size){
          for(j in 1:mediator_size){
            R[i,j]=a[i]*a[j]
          }
        }
        R <- cbind(a,R)
        R <- rbind(c(0,a),R)
        diag(R)=diag(R)+1
        # generate x and m by equation 4
        mu <- rep(0, mediator_size+1)
        X <- mvrnorm(sample_size, mu, R)
        beta = c(r, b)
        beta <- matrix(beta, ncol = 1)
        pr = 1/(1+exp(-(alpha1 + X %*% beta)))  # generate y by equation 2
        y = rbinom(n = sample_size, size = 1, prob = pr)
        # generate noise mediators
        noise = matrix(rnorm(sample_size*(p - mediator_size)),nrow=sample_size,ncol=p - mediator_size)
        M=cbind(X[,2:(mediator_size+1)],noise)
        MT=X[,2:(mediator_size+1)]
        colnames(M)=paste0('V', 1:ncol(M))
        X=X[,1]
        dat <- data.frame(y, X, M)
        return(list(dat=dat,MT=MT, pr=pr))
      }}
    else{
      x = rnorm(n = sample_size)           # generate x
      #specify covariance matrix R among x and m
      R <- matrix(0, ncol = mediator_size,nrow = mediator_size)
      for(i in 1:mediator_size){
        for(j in 1:mediator_size){
          R[i,j]=a[i]*a[j]
        }
      }
      R <- cbind(a,R)
      R <- rbind(c(0,a),R)
      diag(R)=diag(R)+1
      # generate x and m by equation 4
      mu <- rep(0, mediator_size+1)
      m1=matrix(rnorm(sample_size*nm1), nrow=sample_size, ncol=nm1)
      X <- scale(cbind(mvrnorm(sample_size, mu, R),m1))
      beta = c(r, b)
      beta <- matrix(beta, ncol = 1)
      pr = 1/(1+exp(-(alpha1 + X %*% beta)))  # generate y by equation 2
      y = rbinom(n = sample_size, size = 1, prob = pr)
      # generate noise mediators
      noise = matrix(rnorm(sample_size*(p - mediator_size-nm1)),nrow=sample_size,ncol=p - mediator_size-nm1)
      M=cbind(X[,2:(mediator_size+nm1+1)],noise)
      MT=X[,2:(mediator_size+1)]
      colnames(M)=paste0('V', 1:ncol(M))
      X=X[,1]
      dat <- data.frame(y, X, M)
      return(list(dat=dat,MT=MT, pr=pr))
    }
  }}

