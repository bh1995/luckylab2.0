#' Linear Regression.
#' 
#' @param formula A formula.
#' @param data A data frame.
#' @return An object of class linreg.
#' @export
#' 
linreg <- function(formula, data){
  X <- model.matrix(formula,data)
  y <- data[,all.vars(formula)[1]]
  n <- length(y)
  #QR decomposition
  QR <<- function(X){
    # Empty U matrix
    U <<- matrix(0,dim(X)[1],dim(X)[2])
    
    # Fill U matrix
    for (i in 1:(dim(X)[2])){
      U[,i] <- cbind(X[,i]) 
      if (dim(X)[2]>1){
        for (j in 1:(dim(X)[2]-1)){
          if (j+1>i) {break}
          U[,i] <- U[,i] - as.numeric(t(U[,j])%*%X[,i]/t(U[,j])%*%U[,j])*U[,j]
        }
      }
    }
    
    # Calculate the e-values and put them in the Q matrix.
    Q <<- apply(U, 2, function(x)(x/sqrt((sum(x^2)))))
    
    # Empty R matrix
    R <<- matrix(0,dim(X)[2],dim(X)[2])
    
    # Fill R matrix
    for (i in 1:(dim(X)[2])){
      # Start j from (dim(X)[2]) so that lower triangular of R = 0
      for (j in (dim(X)[2]):1){
        if (j<i) {break}
        R[i,j] <- t(Q[,i])%*%X[,j]
      }
    }
    return(list(Q,R))
  }
  
  Q <- QR(X)[[1]]  
  R <- QR(X)[[2]]

  multreg <- function(Q, R, y){
      b <- (solve(R)) %*% t(Q) %*% y
      fitted <- Q%*%R%*%b
      resid <- y - fitted
      df <- nrow(Q)-ncol(Q)
      residvar <- t(resid)%*%resid/df
      varb <- solve(t(R)%*%R) * as.numeric(residvar)
      tval <- b/sqrt(diag(varb))
      pval <- 1-pt(tval,df)
      return(list(b, fitted, resid, df, residvar, varb, tval, pval))
  }
  
  output <<- multreg(Q, R, y)
  

  coeff <<- c(output[[1]])
  names(coeff) <<- colnames(X)
  coeff <- c(output[[1]])
  names(coeff) <- colnames(X)
  
  dataname <- deparse(substitute(data)) #for the print methods
  reg <<- list(formula, coeff, output[[2]], output[[3]], output[[6]], output[[7]], output[[8]], output[[4]], dataname, n)
  names(reg) <<- c("formula","coefficients","fitted","residuals","varcoef","t-values","p-values","df", "dataname", "n")
  class(reg) <<- "linreg" 
  return(reg)
}

#' Print
#'
#' \code{print.linreg} Prints out the coefficients and coefficient names.
#'
#' This function returns the coefficients and coefficient names stored in 
#' the object linreg of S3 class.
#'
#' @export
#' @param x A class object of Linear Regression
#' @return The coefficients and coefficient names.
print <- function(x){UseMethod("print",x)}
#' @export
print.linreg <- function(x){
  a <- as.character(x$dataname)
  b <- as.character(x$formula)
  
  # cat("Call:\n", "linreg(", c )
  # 
  # cat(")",",","data =", a, "\n","Coefficients:\n")
  
  cat("Call:")
  cat("\n")
  # formula_print<- paste0("linreg(formula = ","",b,","," data = ","",a ,")","\n","\n","Coefficients:\n", sep=" " )
  # cat(formula_print)
  # round(x$coefficients, digits=3)
  # x$coefficients
}
# ***test <- linreg(formula, iris) <- must enter "iris" instead of "data"***



#' Residuals
#'
#' \code{resid.linreg} Returns the vector of residuals e
#'
#' This function returns the vector of residuals e stored in 
#' the object linreg of S3 class.
#'
#' @export
#' @param x An object of class linreg.
#' @return The vector of residuals e.
resid <- function(x){UseMethod("resid",x)}
#' @export
resid.linreg <- function(x,...){
  return(x$residuals)
}


#' Fitted values
#'
#' \code{pred.linreg} Returns the predicted values.
#'
#' This function returns the predicted values stored in 
#' the object linreg of S3 class.
#'
#' @export
#' @param x An object of class linreg.
#' @return Predicted values yhat.
pred <- function(x){UseMethod("pred",x)}
#' @export
pred.linreg <- function(x,...){
  return(x$fitted)
}


#' Regression coefficients
#'
#' \code{coef.linreg} Returns the coefficients as a named vector.
#'
#' This function returns the coefficients as a named vector stored in 
#' the object linreg of S3 class.
#'
#' @export
#' @param x An object of class linreg.
#' @return Coefficients as a named vector.
coef <- function(x){UseMethod("coef",x)}
#' @export
coef.linreg <- function(x,...){
  return(x$coefficients)
}


#' Regression summary
#'
#' \code{summary.linreg} Outouts a summary of the calculated results from the regression model
#'
#' @export
#' @param x Which is an object of class linreg
#' @return summary of calculated values
summary <- function(x){UseMethod("summary",x)}

#' @export
summary.linreg <- function(x){
  printSum <- function(b, bvar, tval, pval, n){
    cat(sprintf("Coefficients\t\tSdError\t\tTvalue\t\tPvalue\n"))
  for (i in 1:n){
    
    cat(sprintf("%-15s\t\t%6.6f\t%6.6f\t%6.6f\n",strtrim(names(b[i]),15),bvar[i],tval[i],pval[i]))
  }
  cat(sprintf("Estimated error variance: %f\n",var(x$resid)))
  cat(sprintf("Degrees of freedom: %d\n",x$df))
  
}
  n <- length(x$coefficients)
  coeff <- x$coefficients
  bvar <- x$bvar
  T <- x$tval
  P <- x$pval
  printSummary(coefficients,bvar,T,pval,n)
}

