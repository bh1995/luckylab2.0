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
      bvar <- c(residvar)*diag(solve(t(X)%*%X))
      return(list(b, fitted, resid, df, residvar, varb, t, pval,bvar))
  }
  
  output <- multreg(Q, R, y)
  

  coeff <<- c(output[[1]])
  names(coeff) <<- colnames(X)
  coeff <- c(output[[1]])
  names(coeff) <- colnames(X)
  
  dataname <- deparse(substitute(data)) #for the print methods
  reg <<- list(formula, coeff, output[[2]], output[[3]], output[[6]], output[[7]], output[[8]], output[[4]], dataname, X, output[[9]])
  names(reg) <<- c("formula","coefficients","fitted","residuals","varcoef","t-values","p-values","df", "dataname", "X","bvar")
  class(reg) <<- "linreg" 
  return(reg)
}

#' Regression summary
#'
#' \code{summary.linreg} Outputs a summary of the calculated results from the regression model.
#'
#' @export
#' @param x An object of class linreg.
#' @return Summary of calculated values.
summary <- function(x){UseMethod("summary",x)}
#' @export
summary.linreg <- function(x){
  a <- as.character(x$dataname)
  format_print <- format(x$formula)
  cat("Call:\n", "linreg(formula = ", format_print ,","," data = ","",a,")","\n","\n", sep="")
  
  cat("Coefficients", "Standard Error", "T values", "P values\n")
  n <- length(x$coefficients)
  for( i in 1:n){
    cat(names(x$coefficients)[i], round(x$coefficients[i], digits =3) , 
                                  round(x$bvar[i], digits = 3),
                                  round(x$'t-values'[i], digits=3),
                                  round(x$'p-values'[i], digits =3), "***","\n")
  }
  cat("\n")
  res <- as.vector(x$residuals)
  dff <- as.vector(x$df)
  cat("Residual standard error: ", var(res)," on ", dff, " degrees of freedom", sep="" )
}



#' Print
#'
#' \code{print.linreg} Prints out the coefficients and coefficient names.
#'
#' This function returns the coefficients and coefficient names stored in 
#' the object linreg of S3 class.
#'
#' @export
#' @param x An object of class linreg.
#' @return The coefficients and coefficient names.
print <- function(x){UseMethod("print",x)}
#' @export
print.linreg <- function(x){
  a <- as.character(x$dataname)
  format_print <- format(x$formula)
 
  cat("Call:\n", "linreg(formula = ", format_print ,","," data = ","",a,")","\n","\n","Coefficients:\n", sep="")
  print(x$coefficient)
  
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

#' Plots
#'
#' \code{plot.linreg} Prints residual plots.
#'
#' This function plots residual plots.
#'
#' @export
#' @param x An object of class linreg.
#' @return Residual plots.
plot <- function(x){UseMethod("plot",x)}
#' @export
#' @import ggplot2 gridExtra
#' 
plot.linreg <- function(x){
  # First plot
  residfit <- function(formula, fitted, resid){
    a <- as.character(as.expression(formula))
    out <- tail(order(abs(resid)),3)
    plot <- ggplot(data=NULL,aes(x=fitted, y=resid)) +
      geom_point()  +
      geom_text(aes(fitted[out], resid[out]), label = out, size=3, hjust = 1.2) +
      stat_summary(fun.y="mean", geom="line", aes(fitted), color = "red") +
      labs(y = "Residuals", x = paste0("Fitted values \nlinreg(",a,")")) +
      ggtitle("Residuals vs Fitted") +
      theme(plot.title = element_text(hjust = 0.5))
    return(plot)
  }
  
  # Second plot
  scaleloc <- function(formula, fitted, resid){
    standardized <- sqrt(abs(resid-mean(resid))/sd(resid))
    a <- as.character(as.expression(formula))
    out <- tail(order(abs(resid)),3)
    plot <- ggplot(data=NULL,aes(x=fitted, y=standardized)) +
      geom_point()  +
      geom_text(aes(fitted[out], standardized[out]), label = out, size=3, hjust = 1.2) +
      stat_summary(fun.y="mean", geom="line", aes(fitted), color = "red") +
      labs(y = expression(sqrt(abs("Standardized residuals"))), x = paste0("Fitted values \nlinreg(",a,")")) +
      ggtitle("Scale-Location") +
      theme(plot.title = element_text(hjust = 0.5))
    return(plot)
  }
  # Combine two graphs into one plot
  grid.arrange(residfit(formula=x$formula, fitted=x$fitted, resid=x$residuals), 
               scaleloc(formula=x$formula, fitted=x$fitted, resid=x$residuals), ncol = 2)
}













