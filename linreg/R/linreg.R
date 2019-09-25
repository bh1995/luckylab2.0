#' Linear Regression.
#' 
#' @param formula A object of class "formula".
#' @param data A data frame.
#' @return An S3 object of class linreg.
#' @export
#' 
linreg <- function(formula, data){
  X <- model.matrix(formula,data)
  y <- data[,all.vars(formula)[1]]
  
  #QR decomposition
  QR <- function(X){
    # Empty U matrix
    U <- matrix(0,dim(X)[1],dim(X)[2])
    
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
    Q <- apply(U, 2, function(x)(x/sqrt((sum(x^2)))))
    
    # Empty R matrix
    R <- matrix(0,dim(X)[2],dim(X)[2])
    
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
  
  output <- multreg(Q, R, y)
  
  coeff <- c(output[[1]])
  names(coeff) <- colnames(X)
  
  dataname <- deparse(substitute(data))
  
  reg <- list(formula, coeff, output[[2]], output[[3]], output[[6]], output[[7]], output[[8]], output[[4]], dataname)
  names(reg) <- c("formula","coefficients","fitted","residuals","varcoef","t-values","p-values","df", "dataname")
  class(reg) <- "linreg" 
  return(reg)
}

#' Print
#'
#' \code{print.linreg} print out the coefficients and coefficient names
#'
#' This function returns the coefficients and coefficient names stored in 
#' the object linreg of S3 class.
#'
#' @export
#' @param x A class object of Linear Regression
#' @return the coefficients and coefficient names
print <- function(x){UseMethod("print",x)}
#' @export
print.linreg <- function(x,...){
  a <- as.character(x$dataname)
  b <- as.character(x$formula)
  cat("Call:\n", "linreg(",b[2],b[1],b[3], ")",",","data =", a, "\n\n","Coefficients:\n" )
  x$coefficients
}

# test <- linreg(formula, iris) <- must enter "iris" instead of "data"



#' Residuals
#'
#' \code{resid.linreg} return the vector of residuals e
#'
#' This function returns the vector of residuals e stored in 
#' the object linreg of S3 class.
#'
#' @export
#' @param x A class object of Linear Regression
#' @return the vector of residuals e
resid <- function(x){UseMethod("resid",x)}
#' @export
resid.linreg <- function(x,...){
  return(x$residuals)
}


#' Fitted values
#'
#' \code{pred.linreg} Returns the predicted values yhat
#'
#' This function returns  the predicted values y stored in 
#' the object linreg of S3 class.
#'
#' @export
#' @param x  A class object of Linear Regression
#' @return A predicted values yhat.
pred <- function(x){UseMethod("pred",x)}
#' @export
pred.linreg <- function(x,...){
  return(x$fitted)
}


#' Regression coefficients
#'
#' \code{coef.linreg} Returns the coefficients as a named vector.
#'
#' This function returns the the coefficients as a named vector stored in 
#' the object linreg of S3 class.
#'
#' @export
#' @param x A class object of Linear Regression
#' @return A coefficients as a named vector
coef <- function(x){UseMethod("coef",x)}
#' @export
coef.linreg <- function(x,...){
  return(x$coefficients)
}


