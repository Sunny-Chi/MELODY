#' Partial R-Squared Simulation for Logistic Regression
#'
#' Conducts simulations to assess the partial R-squared values from logistic regression analyses, focusing on the contribution of selected covariates to the model's explanatory power.
#'
#' @param Q Number of simulations to perform.
#' @param sample_size Number of samples per simulation.
#' @param nsis Parameter for Stepwise Information Criterion in mediator selection.
#' @param mediator_size Number of mediator variables.
#' @param p Total number of predictors.
#' @param a Coefficients for mediator interaction.
#' @param r Coefficient for the primary predictor.
#' @param b Coefficients for the mediators.
#' @param alpha1 Intercept for the logistic model.
#' @param Rsmethod Method for calculating R-squared (e.g., "McFadden").
#' @param nm2 Number of noise mediators.
#' @param trainrate Training set proportion.
#' @param FDR False Discovery Rate for selection.
#' @param select Flag to perform selection.
#' @param gamma Gamma parameter for penalization.
#' @param iter Iteration flag for selection process.
#' @param g Coefficient for additional predictors.
#' @param l Coefficient modeling the effect of Z on M.
#' @return A list of simulation results, including R-squared values and selection metrics.
#' @examples
#' r2_results <- r2partial_sim(Q = 10)
#' @export
r2partial_sim= function(Q=Q,sample_size = sample_size, nsis=nsis,mediator_size = mediator_size, p = p, a = a, r = r, b = b, alpha1 = alpha1,Rsmethod=Rsmethod,nm2=nm2,trainrate=trainrate,FDR=FDR,select=select,gamma=gamma,iter=iter,g=g,l=l){
  pr=rep(NA,Q)
  R2_Z=rep(NA,Q)
  R2_XZ=rep(NA,Q)
  R2_MZ=rep(NA,Q)
  R2_XMZ=rep(NA,Q)
  R2_med  <- rep(NA,Q)
  SOS=rep(NA,Q)
  product= matrix(NA,Q,5)
  colnames(product)=c("ab/c","ab","c-r","(c-r)/c","c")
  FP= matrix(NA,Q)
  FPnm2= matrix(NA,Q)
  TP= matrix(NA,Q)
  for(q in 1:Q){
    print(q)
    set.seed(q)
    df = sim_logistic_data (sample_size = sample_size, mediator_size = mediator_size, p = p, a = a, r = r, b =b,alpha1=alpha1,cov=cov,g=g,l=l)

    #true mediators' name
    if(is.null(nm2)){
      cc=paste0('V', 1:mediator_size)
    }else{
      cc=paste0('V', 1:(mediator_size-nm2))
    }

    Y=df$dat$y
    X=df$dat$X
    Z=df$dat$Z
    M=scale(df$dat[,c(paste0('V', 1:p))])

    if(is.null(select)==FALSE){
      train=1:(sample_size*trainrate)
      # mediator selection
      #step 1&2
      yt=Y[train]
      xt=scale(cbind(X,Z,M)[train,])
      model=SIS(xt, yt, family='binomial', nsis=nsis,penalty = 'MCP', concavity.parameter =gamma,tune='cv',iter=iter)
      #selected by step 1
      #sis.ix0 is index, exposure is 1 cov is 2 in sis.ix0 so we minus 2 for index to get the mediators' index
      SIS_M=model$sis.ix0-2
      SIS_M=SIS_M[SIS_M>0]
      SIS_M=colnames(M)[SIS_M]
      print(SIS_M)
      #selected by step 2
      MCP_M<-model$ix-2
      MCP_M=MCP_M[MCP_M>0]
      MCP_M<-colnames(M)[MCP_M]
      M_MCP <- M[train, MCP_M]
      print(MCP_M)
      if(sum(model$ix>2)==0) {
        R2_med[q]  <- 0
        SOS[q]=0
        product[q,]=c(0,0,NA,0,0)
        coef_a=0
        coef_b=0
        next
      }


      #step 3
      if(is.null(FDR)==FALSE){

        cal_alpha_simple<-function(x){
          data2<-data.frame(Med=x,envir=X[train])
          l<-summary(stats::lm('Med ~.', data=data2))
          invalphap<-(l$coefficients)['envir','Pr(>|t|)']
          return(invalphap)
        }
        # keep the top candidates for ultra high dimensional data
        if(length(MCP_M)==1) {inv.alpha.p<-cal_alpha_simple(M_MCP)}else{
          inv.alpha.p<-apply(M_MCP,2, cal_alpha_simple)}
        #order <- order(stats::p.adjust(inv.alpha.p, method='fdr'))[1:round(init.cutoff *ncol(Med)) ]
        FDR_M=names(which(stats::p.adjust(inv.alpha.p, method='fdr')<FDR))
        print(FDR_M)
        M=scale(M[-train,FDR_M])
        pp=sum(cc %in% unlist(FDR_M))
        FP[q]=length(FDR_M)-pp
        TP[q]=pp
        FPnm2[q]=sum(paste0('V', 1:mediator_size)%in% unlist(FDR_M))-pp #number of false selected m2
        if(length(FDR_M)==0) {
          R2_med[q]  <- 0
          SOS[q]=0
          product[q,]=c(0,0,NA,0,0)
          coef_a=0
          coef_b=0
          next
        }
      }else{
        M=scale(M[-train,MCP_M])
        pp=sum(cc %in% unlist(MCP_M))
        FP[q]=length(MCP_M)-pp
        TP[q]=pp
        FPnm2[q]=sum(paste0('V', 1:mediator_size)%in% unlist(MCP_M))-pp #number of false selected m2
      }

      Y=Y[-train]
      X=X[-train]
      Z=Z[-train]
    }


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

    pr[q]=mean(df$pr)
    R2_med[q]  <- (R2_XZ[q]  + R2_MZ[q]  - R2_XMZ[q] - R2_Z[q])/(1-R2_Z[q])
    SOS[q]=R2_med[q]/((R2_XZ[q]-R2_Z[q])/(1-R2_Z[q]))
    product[q,]=c(coef_a%*%coef_b/coef_c,coef_a%*%coef_b,coef_c-coef_r,(coef_c-coef_r)/coef_c,coef_c)
  }
  if(is.null(nm2)){
    tpr=TP/mediator_size
    fpr=FP/(p-mediator_size)
  }else{
    tpr=TP/(mediator_size-nm2)
    fpr=FP/(p-mediator_size+nm2)
  }
  return(list(pr=pr,R2_Z=R2_Z,R2_XZ=R2_XZ,R2_MZ=R2_MZ,R2_XMZ=R2_XMZ,R2_med=R2_med,SOS=SOS,product=product,fpr=fpr,tpr=tpr,fp=FP,tp=TP,fpnm2=FPnm2,coef_a=coef_a,coef_b=coef_b))
}
