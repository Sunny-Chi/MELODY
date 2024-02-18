#' Mediation model and estimation function
#'
#' This function simulates logistic regression data and analyzes it using various models
#' to calculate R-squared values and other statistics. It can handle data with or without
#' a predefined covariance structure and allows for customization of simulation parameters.
#'
#' @param Q Number of simulations to run.
#' @param sample_size Number of samples per simulation.
#' @param mediator_size Number of mediator variables.
#' @param p Total number of predictors.
#' @param a Coefficients for mediators.
#' @param r Coefficient for the primary predictor.
#' @param b Coefficients for the mediators.
#' @param alpha1 Intercept.
#' @param Rsmethod Method for calculating R-squared (e.g., "Tjur").
#' @param cov Optional covariance matrix.
#' @param g Additional parameter for `sim_logistic_cov`.
#' @param l Additional parameter for `sim_logistic_cov`.
#' @return A list containing simulation results, including R-squared values, coefficients, and probabilities.
#' @examples
#' result <- sim_fun()
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom glmnet glm
#' @importFrom DescTools PseudoR2

sim_fun = function(Q=500, sample_size = 2000, mediator_size = 2, p = 1000, a = rep(0,2), r = 1, b = rep(0,2), alpha1 = 0,Rsmethod="Tjur",cov=NULL,g=g,l=l) {
  if(is.null(cov)==FALSE){
    pr=rep(NA,Q)
    R2_Z=rep(NA,Q)
    R2_XZ=rep(NA,Q)
    R2_MZ=rep(NA,Q)
    R2_XMZ=rep(NA,Q)
    R2_med  <- rep(NA,Q)
    SOS=rep(NA,Q)
    product= matrix(NA,Q,5)
    colnames(product)=c("ab/c","ab","c","c-r","(c-r)/c")
    for(q in 1:Q){
      print(q)
      set.seed(q)
      df = sim_logistic_data (sample_size = sample_size, mediator_size = mediator_size, p = p, a = a, r = r, b =b,alpha1=alpha1,cov=cov,g=g,l=l)
      Y=df$dat$y
      X=scale(df$dat$X)
      Z=df$dat$Z
      #M=scale(df$dat[,c(paste0('V', 1:p))])
      M=scale(df$MT)
      glm_XMZ <- glm(Y ~ Z+X+M, family=binomial)
      R2_XMZ[q]=PseudoR2(glm_XMZ,Rsmethod)
      coef_b = glm_XMZ$coefficients[-c(1,2,3)]
      if(max(coef_b)>50) {
        R2_med[q]  <- 0
        SOS[q]=0
        product[q,]=c(0,0,NA,0,0)
        coef_a=0
        coef_b=0
        next
      }
      coef_r = glm_XMZ$coefficients[3]
      glm_MZ <- glm(Y ~ Z+M, family=binomial)
      R2_MZ[q]=PseudoR2(glm_MZ,Rsmethod)
      med1 <- lm(M ~ Z+X)
      coef_a = med1$coefficients
      if(is.null(dim(coef_a))){
        coef_a = coef_a[3]
      }else{coef_a = coef_a[3,]}
      glm_XZ <- glm(Y ~ Z+X, family=binomial)
      coef_c = glm_XZ$coefficients[3]
      R2_XZ[q]=PseudoR2(glm_XZ,Rsmethod)
      glm_Z <- glm(Y ~ Z, family=binomial)
      R2_Z[q]=PseudoR2(glm_Z,Rsmethod)
      R2_med[q]  <- (R2_XZ[q]  + R2_MZ[q]  - R2_XMZ[q] - R2_Z[q])/(1-R2_Z[q])
      SOS[q]=R2_med[q]/((R2_XZ[q]-R2_Z[q])/(1-R2_Z[q]))
      product[q,]=c(coef_a%*%coef_b/coef_c,coef_a%*%coef_b,coef_c,coef_c-coef_r,(coef_c-coef_r)/coef_c)
      pr[q]=mean(df$pr)
    }
    return(list(pr=pr,R2_Z=R2_Z,R2_XZ=R2_XZ,R2_MZ=R2_MZ,R2_XMZ=R2_XMZ,R2_med=R2_med,SOS=SOS,product=product,coef_a=coef_a,coef_b=coef_b))
  }else{
    pr=rep(NA,Q)
    R2_X=rep(NA,Q)
    R2_M=rep(NA,Q)
    R2_XM=rep(NA,Q)
    R2_med  <- rep(NA,Q)
    SOS=rep(NA,Q)
    product= matrix(NA,Q,5)
    colnames(product)=c("ab/c","ab","c","c-r","(c-r)/c")

    for(q in 1:Q){
      print(q)
      set.seed(q)
      df = sim_logistic_data (sample_size = sample_size, mediator_size = mediator_size, p = p, a = a, r = r, b =b,alpha1=alpha1,cov=cov)
      Y=df$dat$y
      X=scale(df$dat$X)
      #M=scale(df$dat[,c(paste0('V', 1:p))])
      M=scale(df$MT)

      glm_XM <- glm(Y ~ X+M, family=binomial)
      R2_XM[q]=PseudoR2(glm_XM,Rsmethod)
      coef_b = glm_XM$coefficients[-c(1,2)]
      if(max(coef_b)>50) {
        R2_med[q]  <- 0
        SOS[q]=0
        product[q,]=c(0,0,NA,0,0)
        coef_a=0
        coef_b=0
        next
      }
      coef_r = glm_XM$coefficients[2]

      glm_M <- glm(Y ~ M, family=binomial)
      R2_M[q]=PseudoR2(glm_M,Rsmethod)

      med1 <- lm(M ~ X)
      coef_a = med1$coefficients[2,]

      glm_X <- glm(Y ~ X, family=binomial)
      coef_c = glm_X$coefficients[2]
      R2_X[q]=PseudoR2(glm_X,Rsmethod)

      pr[q]=mean(df$pr)
      R2_med[q]  <- R2_X[q]  + R2_M[q]  - R2_XM[q]
      SOS[q]=R2_med[q]/R2_X[q]
      product[q,]=c(coef_a%*%coef_b/coef_c,coef_a%*%coef_b,coef_c,coef_c-coef_r,(coef_c-coef_r)/coef_c)
    }
    return(list(pr=pr,R2_X=R2_X,R2_M=R2_M,R2_XM=R2_XM,R2_med=R2_med,SOS=SOS,product=product,coef_a=coef_a,coef_b=coef_b))
  }}
