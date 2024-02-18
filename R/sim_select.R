#' Mediator selection function
#'
#' This function conducts selective simulations to evaluate the impact of various mediators and predictors on logistic regression outcomes. It allows for mediator selection based on specified criteria and calculates performance metrics.
#'
#' @param Q Number of simulations.
#' @param sample_size Number of observations per simulation.
#' @param nsis Number of selections in Stepwise Information Criterion.
#' @param mediator_size Number of mediators.
#' @param p Total number of predictors.
#' @param a Coefficients for mediator interaction.
#' @param r Coefficient for the primary predictor.
#' @param b Coefficients for the mediators.
#' @param alpha1 Intercept for the logistic model.
#' @param Rsmethod Method for R-squared calculation.
#' @param nm2 Number of noise mediators to include.
#' @param nm1 Number of additional noise mediators.
#' @param trainrate Proportion of data used for training.
#' @param cov Optional covariance matrix for data simulation.
#' @param FDR False Discovery Rate for mediator selection.
#' @param select Boolean to activate selection process.
#' @param gamma Gamma parameter for penalization in mediator selection.
#' @param iter Boolean to iterate selection process.
#' @param g Additional parameter for `sim_logistic_cov`.
#' @param l Additional parameter for `sim_logistic_cov`.
#' @return A list of simulation results including R-squared values, selection performance metrics, and coefficients.
#' @examples
#' results <- sim_select()
#' @export

sim_select = function(Q=500, sample_size = 2000, nsis=2000/(4*log(2000)),mediator_size = 2, p = 1000, a = rep(0,2), r = 1, b = rep(0,2), alpha1 = 0,Rsmethod="McFadden",nm2=NULL,nm1=NULL,trainrate=0.75,cov=NULL,FDR=0.2,select=TRUE,gamma=3,iter=FALSE,g=NULL,l=NULL) {
  if(is.null(cov)==FALSE){
    return(r2partial_sim(Q=Q,sample_size = sample_size, nsis=nsis,mediator_size = mediator_size, p = p, a = a, r = r, b = b, alpha1 = alpha1,Rsmethod=Rsmethod,nm2=nm2,trainrate=trainrate,FDR=FDR,select=select,gamma=gamma,iter=iter,g=g,l=l))
  }else{

    pr=rep(NA,Q)
    R2_X=rep(NA,Q)
    R2_M=rep(NA,Q)
    R2_XM=rep(NA,Q)
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
      df = sim_logistic_data (sample_size = sample_size, mediator_size = mediator_size, p = p, a = a, r = r, b =b,alpha1=alpha1,cov=cov,nm1=nm1)

      #true mediators' name
      if(is.null(nm2)){
        cc=paste0('V', 1:mediator_size)
      }else{
        cc=paste0('V', 1:(mediator_size-nm2))
      }



      Y=df$dat$y
      X=scale(df$dat$X)
      M=scale(df$dat[,c(paste0('V', 1:p))])

      if(is.null(select)==FALSE){
        train=1:(sample_size*trainrate)
        # mediator selection
        #step 1&2
        yt=Y[train]
        xt=scale(cbind(X,M)[train,])
        model=SIS(xt, yt, family='binomial', nsis=nsis,penalty = 'MCP', concavity.parameter =gamma,tune='cv',iter=iter)
        #selected by step 1
        #sis.ix0 is index, exposure is 1 in sis.ix0 so we minus 1 for index to get the mediators' index

        SIS_M=colnames(M[train, model$sis.ix0-1])
        print(SIS_M)
        #selected by step 2
        MCP_M<-colnames(M)[model$ix-1]
        M_MCP <- M[train, MCP_M]
        print(MCP_M)
        if(sum(model$ix>1)==0) {
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
        X=scale(X[-train])
      }

      glm_XM <- glm(Y ~ X+M, family=binomial)
      R2_XM[q]=PseudoR2(glm_XM,Rsmethod)
      coef_b = glm_XM$coefficients[-c(1,2)]

      if(max(coef_b,na.rm=TRUE)>50) {
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
      coef_a = med1$coefficients
      if(is.null(dim(coef_a))){
        coef_a = coef_a[2]
      }else{coef_a = coef_a[2,]}

      glm_X <- glm(Y ~ X, family=binomial)
      coef_c = glm_X$coefficients[2]
      R2_X[q]=PseudoR2(glm_X,Rsmethod)

      pr[q]=mean(df$pr)
      R2_med[q]  <- R2_X[q]  + R2_M[q]  - R2_XM[q]
      SOS[q]=R2_med[q]/R2_X[q]
      product[q,]=c(coef_a%*%coef_b/coef_c,coef_a%*%coef_b,coef_c-coef_r,(coef_c-coef_r)/coef_c,coef_c)
    }

    if(is.null(nm2)){
      tpr=TP/mediator_size
      fpr=FP/(p-mediator_size)
    }else{
      tpr=TP/(mediator_size-nm2)
      fpr=FP/(p-mediator_size+nm2)
    }
    return(list(pr=pr,R2_X=R2_X,R2_M=R2_M,R2_XM=R2_XM,R2_med=R2_med,SOS=SOS,product=product,fpr=fpr,tpr=tpr,fp=FP,tp=TP,fpnm2=FPnm2,coef_a=coef_a,coef_b=coef_b))
  }}
