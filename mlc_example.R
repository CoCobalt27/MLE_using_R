# Basic ####
# Packages
library(tidyr)
library(dplyr)
library(magrittr)
library(olsrr)

# Data import ####
load("./bc_example.RData")

# Positively constrained multiple linear regression (pMLR) ####
pMLR<-function(onesite, method, intercept){
  
  if(!method %in% c("OLS", "MLE")) stop("Please input a proper method name: OLS or MLE!")
  if(!is.logical(intercept)) stop("Please use a logical variable in intercept item.")
  
  # OLS w/o intercept
  func1<-function(k){
    EC_pred<-onesite$resid*k[1]+onesite$trans*k[2]+onesite$power*k[3]+onesite$indus*k[4]+
             onesite$other*k[5]+onesite$burn*k[6]+onesite$dust*k[7]
    sum((onesite$bc_obs-EC_pred)^2)
  }
  
  # OLS w/ intercept
  func2<-function(k){
    EC_pred<-onesite$resid*k[1]+onesite$trans*k[2]+onesite$power*k[3]+onesite$indus*k[4]+
             onesite$other*k[5]+onesite$burn*k[6]+onesite$dust*k[7]+k[8]
    sum((onesite$bc_obs-EC_pred)^2)
  }
  
  # MLE w/o intercept
  func3<-function(k){
    EC_pred<-onesite$resid*k[1]+onesite$trans*k[2]+onesite$power*k[3]+onesite$indus*k[4]+
             onesite$other*k[5]+onesite$burn*k[6]+onesite$dust*k[7]
    sum((log(onesite$bc_obs)-log(EC_pred))^2)
  }
  
  # MLE w/ intercept
  func4<-function(k){
    EC_pred<-onesite$resid*k[1]+onesite$trans*k[2]+onesite$power*k[3]+onesite$indus*k[4]+
             onesite$other*k[5]+onesite$burn*k[6]+onesite$dust*k[7]+k[8]
    sum((log(onesite$bc_obs)-log(EC_pred))^2)
  }
  
  
  # Parameter estimation
  if(intercept){
    if(method == "OLS"){
      optim_result<-optim(rep(0.5, 8), func2, method = "L-BFGS-B", lower = rep(0.000001, 8))
    }else{
      optim_result<-optim(rep(0.5, 8), func4, method = "L-BFGS-B", lower = rep(0.000001, 8))
    }
    optim_result$par->para
    optim_result$value->Q
  }else{
    if(method == "OLS"){
      optim_result<-optim(rep(0.5, 7), func1, method = "L-BFGS-B", lower = rep(0.000001, 7))
    }else{
      optim_result<-optim(rep(0.5, 7), func3, method = "L-BFGS-B", lower = rep(0.000001, 7))
    }
    c(optim_result$par, 0)->para
    optim_result$value->Q
  }
  names(para)<-names(onesite)[3:10]
  
  # Fitted value
  onesite$fit_value<-onesite$resid*para[1]+onesite$trans*para[2]+onesite$power*para[3]+
                     onesite$indus*para[4]+onesite$other*para[5]+onesite$burn*para[6]+
                     onesite$dust*para[7]+para[8]
  
  # Residuals
  if(method == "OLS"){
    onesite$res<-onesite$bc_obs - onesite$fit_value
  }else{
    onesite$res<-log(onesite$bc_obs / onesite$fit_value)
  }
  
  # Normality test for residuals
  normal_p_value<-ks.test(onesite$res, "pnorm", 0, sd(onesite$res), exact = FALSE)$p.value
  
  # Homoscedasticity test
  onesite%>%
    mutate(g = res^2/sum(res^2)*nrow(.))%>%
    lm(g ~ fit_value, .)->lm2
  onesite%>%
    mutate(g = res^2/sum(res^2)*nrow(.))%$%
    sum((g-1)^2)-sum(lm2$residuals^2)->SS
  pchisq(SS/2, 1, lower.tail = FALSE)->homosce_p_value
  
  # PCC, MFE, MFB
  cor(onesite$fit_value, onesite$bc_obs)->PCC
  2 * mean((onesite$fit_value - onesite$bc_obs)/(onesite$fit_value + onesite$bc_obs))->MFB
  2 * mean(abs(onesite$fit_value - onesite$bc_obs)/(onesite$fit_value + onesite$bc_obs))->MFE
  
  # AIC calculation
  if(intercept){
    if(method == "OLS"){
      AIC<-2 * 9 + nrow(onesite) * log(Q / nrow(onesite))
    }else{
      AIC<-2 * 9 + nrow(onesite) * log(Q / nrow(onesite)) + sum(log(onesite$bc_obs))
    }
  }else{
    if(method == "OLS"){
      AIC<-2 * 8 + nrow(onesite) * log(Q / nrow(onesite))
    }else{
      AIC<-2 * 8 + nrow(onesite) * log(Q / nrow(onesite)) + sum(log(onesite$bc_obs))
    }
  }
  
  # Results arranged as a
  list("meta-info" = c("site" = onesite$site[1], "method" = method, "w/ intercept" = intercept),
       "coefficient" = round(para, 3),
       "prediction metrics" = round(c("PPC" = PCC, "MFE" = MFE, "MFB" = MFB),3),
       "residual analysis" = c("K-S test p value" = normal_p_value,
                               "B-P test p value" = homosce_p_value),
       "AIC" = AIC)%>%
    return()
}

pMLR(BC, "MLE", TRUE)