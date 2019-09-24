#' Regression
#' 
#' @param Q is matrix formed via QR decomposition
#' @param R is a triangular matrix formed via QR decomposition
#' @param y is the variable that is to be regressed on 
#' @return list of the calculated coefficients
#' @examples
#' 
#' @export
multreg <- function(Q, R, y){
  b <- (solve(R)) %*% t(Q) %*% y
  fitted <- Q%*%R%*%b
  resid <- y - fitted
  df <- nrow(Q%*%R)-ncol(Q%*%R)
  residvar <- t(resid)%*%resid/df
  varb <- solve(t(R)%*%R) * as.numeric(residvar)
  t <- b/sqrt(diag(varb))
  pval <- 1-pt(t,df)
  return(list(b, fitted, resid, df, residvar, varb, t, pval))
}

