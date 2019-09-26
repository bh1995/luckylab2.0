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
      t <- b/sqrt(diag(varb))
      pval <- 1-pt(t,df)
      return(list(b, fitted, resid, df, residvar, varb, t, pval))
  }
  
  output <<- multreg(Q, R, y)
  
<<<<<<< HEAD
  coeff <<- c(output[[1]])
  names(coeff) <<- colnames(X)
=======
  coeff <- c(output[[1]])
  names(coeff) <- colnames(X)
  
>>>>>>> e0c65301aa25a8bff975b340939a4eb2e24b57b3
  dataname <- deparse(substitute(data)) #for the print methods
  reg <<- list(formula, coeff, output[[2]], output[[3]], output[[6]], output[[7]], output[[8]], output[[4]], dataname, X)
  names(reg) <<- c("formula","coefficients","fitted","residuals","varcoef","t-values","p-values","df", "dataname", "X")
  class(reg) <<- "linreg" 
  return(reg)
}

#' Regression summary
#'
#' \code{plot.linreg} Outouts a summary of the calculated results from the regression model
#'
#' @export
#' @param x Which is an object of class linreg
#' @return summary of calculated values
summary <- function(x){UseMethod("summary",x)}


#' Print
#'
#' \code{print.linreg} Prints out the coefficients and coefficient names.
#'
#' This function returns the coefficients and coefficient names stored in 
#' the object linreg of S3 class.
#'
#' @export
<<<<<<< HEAD
#' @param x A class object of Linear Regression
#' @return the coefficients and coefficient names
print <- function(a){UseMethod("print",a)}
=======
#' @param x An object of class linreg.
#' @return The coefficients and coefficient names.
print <- function(x){UseMethod("print",x)}
>>>>>>> e0c65301aa25a8bff975b340939a4eb2e24b57b3
#' @export

print.linreg <- function(x,...){
<<<<<<< HEAD
  a <<- as.character(x$dataname)
  b <<- as.character(x$formula)
  printformula <- function(x){
    
    return(x$formula)
    
  }
  c <- printformula(x)
  cat("Call:\n", "linreg(", c )
  
  cat( ")",",","data =", a, "\n\n","Coefficients:\n" )
  
  x$coefficients
  
=======
  a <- as.character(x$dataname)
  b <- format(x$formula)
  cat("Call:\n", "linreg(formula = ", b,","," data = ", a, ")","\n\n","Coefficients:\n", sep="" )
  round(x$coefficients, digits=2)
>>>>>>> e0c65301aa25a8bff975b340939a4eb2e24b57b3
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




