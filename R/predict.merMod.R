#' A predict method for merMod models with the bootstrap
#'@description When this package is loaded after loading \code{lme4}, it replaces the predict method for linear mixed-effects models (\code{merMod} objects) with this function. The method provided here adds the argument \code{se.fit} (just like in \code{\link{predict.lm}}). It uses the bootstrap, as implemented in \code{\link{bootMer}} to estimate standard errors and confidence intervals for the predictions. Also useful in conjunction with the \code{\link{visreg}} package, which by default does not plot confidence intervals for \code{merMod} models.
#'@param object An object as returned by \code{\link{lmer}}
#'@param nsim The number of bootstrap replicates. The default is a small number to allow quick testing. A warning is printed when nsim < 100, and it can be set via \code{\link{options}} (see Examples)
#'@param se.fit If TRUE, returns standard error (se.fit) and confidence interval (ci.fit) for the predictions, as components of the returned list.
#'@param alpha Controls the coverage of the confidence interval. 
#'@export
#'@examples 
#'
#'# Fit a linear mixed-effects model 
#'fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#'
#'# Predictions without standard error (fixed effects only)
#'predict(fm1, newdata=data.frame(Days=5), re.form=NA)
#'
#'# Predictions with standard error and confidence interval
#'# Set number of bootstrap replicates first (you should use a larger number)
#'options(nbootsim=20)
#'predict(fm1, newdata=data.frame(Days=5), re.form=NA, se.fit=TRUE)
#'
#'# Add confidence intervals to visreg
#'library(visreg)
#'visreg(fm1, "Days")
#'
#'# Also works with an overlay plot, using visreg
#'# First add an artificial group to the sleepstudy data
#'x <- sort(with(sleepstudy, tapply(Reaction, Subject, mean)))
#'sleepstudy$Group <- as.factor(ifelse(sleepstudy$Subject %in% names(x[1:9]), "A", "B"))
#'fm2 <- lmer(Reaction ~ Days*Group + (Days | Subject), sleepstudy)
#'visreg(fm2, "Days", by="Group", overlay=TRUE)
#'
predict.merMod <- function(object, nsim=getOption("bootnsim"), se.fit=FALSE, alpha=0.05, ...){
  
  if(se.fit){
    if(is.null(nsim))nsim <- 10
    if(nsim < 100)message("Number of bootstrap replicates very low. \n Set to higher value with e.g. options(bootnsim = 500)")
    
    b <- bootMer(object, FUN=function(x)lme4:::predict.merMod(x, ...), nsim=nsim)
    serr <- apply(b$t,2,sd)
    ci <- apply(b$t,2,quantile,probs=c(alpha/2, 1 - alpha/2))
    return(list(fit=b$t0, se.fit=serr, ci.fit=ci))
  } else {
    return(lme4:::predict.merMod(object, ...))
  }
  
}

