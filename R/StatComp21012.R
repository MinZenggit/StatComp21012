#' @useDynLib StatComp21012
#' @importFrom Rcpp sourceCpp
NULL
#> NULL



#' @title Adjusting selection bias and misclassification by IPTW.
#' @description Adjusting selection bias and misclassification by IPTW.
#' @param dt The dataset contains treatment, outcome, covariates, and propensity score.
#' @param Treat Treatment
#' @param Y outcome
#' @param ehat propensity score
#' @param s01 selection parameter P(S=1|Y=0)
#' @param s11 selection parameter P(S=1|Y=1)
#' @param p01 misclassification parameter P(Y*=1|Y=0)
#' @param p11 misclassification parameter P(Y*=1|Y=1)
#' @param strata strata data according to covariates
#' @import stats
#' @return adjusted ATE
#' @examples
#' \dontrun{
#' data(dt)
#' IPTWmethod(dt, Treat, Y, ehat, 1, 1, 0, 1)
#' }
#' @export
IPTWmethod <- function(dt, Treat, Y, ehat, s01, s11, p01, p11, strata = F){
  iptw <- function(dtset){
    nE1 <- mean( dtset$Treat * dtset$Y/ dtset$ehat)
    nE2 <- mean((1-dtset$Treat)*dtset$Y/(1-dtset$ehat))
    A <- (nE1-p01)/s11/(p11-p01)
    B <- (nE2-p01)/s11/(p11-p01)
    E1hat <- A*s01/(1+(s01-s11)*A)
    E2hat <- B*s01/(1+(s01-s11)*B)
    return(E1hat-E2hat)
  }
  if(strata==T){
    dtstrat <- split(x=dt, f=list(dt$X1_ind, dt$X2, dt$X3))
    return(mean(unlist(lapply(dtstrat, iptw))))
  }
  return(iptw(dt))
}


#' @title Adjusting selection bias and misclassification by Double Robust method.
#' @description Adjusting selection bias and misclassification by Double Robust method.
#' @param dt The dataset contains treatment, outcome, covariates, and propensity score.
#' @param object The outcome model
#' @param Treat Treatment
#' @param Y outcome
#' @param ehat propensity score
#' @param s01 selection parameter P(S=1|Y=0)
#' @param s11 selection parameter P(S=1|Y=1)
#' @param p01 misclassification parameter P(Y*=1|Y=0)
#' @param p11 misclassification parameter P(Y*=1|Y=1)
#' @param strata strata data according to covariates
#' @return adjusted ATE
#' @import stats
#' @examples
#' \dontrun{
#' data(dt)
#' DRmethod(dt, lm(Y~X1+X2+X3+Treat, data = dt), Treat, Y, ehat, 1, 1, 0, 1)
#' }
#' @export
DRmethod <- function(dt, object, Treat, Y, ehat, s01, s11, p01, p11, strata = F){
  Drubust <-function(dtset){
    nE1 <- mean( dtset$Treat * (dtset$Y - dtset$ghat1)/ dtset$ehat + dtset$ghat1)
    nE2 <- mean((1-dtset$Treat)*(dtset$Y - dtset$ghat0)/(1-dtset$ehat) + dtset$ghat0)
    A <- (nE1-p01)/s11/(p11-p01)
    B <- (nE2-p01)/s11/(p11-p01)
    E1hat <- A*s01/(1+(s01-s11)*A)
    E2hat <- B*s01/(1+(s01-s11)*B)
    return(E1hat-E2hat)
  }
  dt_treat1 <- dt
  dt_treat1$Treat=1
  dt_treat0 <- dt
  dt_treat0$Treat=0
  dt$ghat1 <- predict.glm(object, newdata = dt_treat1, type = "response")
  dt$ghat0 <- predict.glm(object, newdata = dt_treat0, type = "response")
  if(strata == T){
    dtstrat <- split(x=dt, f=list(dt$X1_ind, dt$X2, dt$X3))
    return(mean(unlist(lapply(dtstrat, Drubust))))
  }
  else
    return(Drubust(dt))
}


#' @title Adjusting selection bias and misclassification by Spline method.
#' @description Adjusting selection bias and misclassification by Spline method.
#' @param dt The dataset contains treatment, outcome, covariates, and propensity score.
#' @param s01 selection parameter P(S=1|Y=0)
#' @param s11 selection parameter P(S=1|Y=1)
#' @param p01 misclassification parameter P(Y*=1|Y=0)
#' @param p11 misclassification parameter P(Y*=1|Y=1)
#' @param df the degree of B-spline
#' @param Boot MCmethod for eastimating the integration, size of the MC sample.
#' @return adjusted ATE
#' @import stats
#' @import splines
#' @examples
#' \dontrun{
#' data(dt)
#' splinemethod(dt, 1, 1, 0, 1)
#' }
#' @export
splinemethod <- function(dt, s01, s11, p01, p11, df = 10, Boot = 1000){
  splinefit <- glm(Y ~ bs(X1, df = df)+ X2 + X3 + Treat,
                   family = binomial(link = "logit"),
                   data = dt)
  ###integrate by Monte Carlo Method
  inte_x1 <- runif(Boot, 0, 80)
  inte_x2 <- rbinom(Boot, 1, 0.5)
  inte_x3 <- sample(c(1:4), Boot, replace = T, prob = rep(0.2, 4))
  indt1 <- data.frame(X1 = inte_x1, X2 = inte_x2, X3 = inte_x3, Treat = 1)
  indt0 <- data.frame(X1 = inte_x1, X2 = inte_x2, X3 = inte_x3, Treat = 0)
  g1 <- predict(splinefit, newdata = indt1, type = "response")
  g0 <- predict(splinefit, newdata = indt0, type = "response")
  A <- (g1-p01)/s11/(p11-p01)
  B <- (g0-p01)/s11/(p11-p01)
  E1 <- A*s01/(1+(s01-s11)*A)
  E0 <- B*s01/(1+(s01-s11)*B)
  return(mean(E1-E0))
}


#' @title Adjusting selection bias and misclassification by Spline method.
#' @description Adjusting selection bias and misclassification by Spline method.
#' @param object The outcome model
#' @param s01 selection parameter P(S=1|Y=0)
#' @param s11 selection parameter P(S=1|Y=1)
#' @param p01 misclassification parameter P(Y*=1|Y=0)
#' @param p11 misclassification parameter P(Y*=1|Y=1)
#' @param Boot MCmethod for eastimating the integration, size of the MC sample.
#' @return adjusted ATE
#' @import stats
#' @examples
#' \dontrun{
#' data(dt)
#' gfit <- glm(formula = Y~X1+X2+X3+Treat, family = binomial(link = "logit"), data = dt)
#' paramethod(gfit, 1, 1, 0, 1)
#' }
#' @export
paramethod <- function(object, s01, s11, p01, p11, Boot = 1000){
  beta_hat <- object$coefficients
  expit <- function(x){exp(x)/(1+exp(x))}
  ### x is a vector of (x1,x2,x3)
  ghat1 <- function(x){
    expit(beta_hat[1] + sum(beta_hat[2:4]*x) + beta_hat[5]*1)
  }
  ghat0 <- function(x){
    expit(beta_hat[1] + sum(beta_hat[2:4]*x) + beta_hat[5]*0)
  }
  A <- function(x){(ghat1(x)-p01)/s11/(p11-p01)}
  B <- function(x){(ghat0(x)-p01)/s11/(p11-p01)}
  Ehat1 <- function(x){A(x)*s01/(1+(s01-s11)*A(x))}
  Ehat0 <- function(x){B(x)*s01/(1+(s01-s11)*B(x))}
  f <- function(x){(Ehat1(x)-Ehat0(x))}
  #integrate by Monte Carlo Method
  inte_x1 <- runif(Boot, 0, 80)
  inte_x2 <- rbinom(Boot, 1, 0.5)
  inte_x3 <- sample(c(1:4), Boot, replace = T, prob = rep(0.25, 4))
  te <- cbind(inte_x1, inte_x2, inte_x3)
  return(mean(apply(te, 1, f)))
}
