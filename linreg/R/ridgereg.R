#'
#' @param formula formula
#' @param data data.frame
#' @param lambda numeric
#' @param x_norm matrix
#' @param beta_ridge numeric
#' @param y_hat numeric
#' @param data_name character
#' @return a ridgereg object.
#' 
#' @importFrom methods new
#' 
#' @export
#' @exportClass ridgereg
#' 
ridgereg <- setRefClass("ridgreg",
                        fields = list(
                          formula = "formula",
                          data = "data.frame",
                          lambda = "numeric",
                          xnorm = "matrix",
                          beta_ridge = "numeric",
                          y_hat = "numeric",
                          data_name = "character"
                        ),
                        methods = list(
                          initialize = function(formula, data, lambda){
                            formula <<- formula
                            data <<- data
                            data_name <<- deparse(substitute(data))
                            lambda <<- lambda
                            
                            # normalization of the covariates
                            
                            x <- model.matrix(object = formula, data = data)
                            for (i in 2:ncol(x)){
                              x[,i] <- ((x[,i]-mean(x[,i]))/ (sqrt(var((x[,i])))))
                            }
                            xnorm <<- x
                            
                            # QR decomposition
                            
                            QR <- qr(x)
                            R <- qr.R(QR)
                            
                            # coeff. beta ridge
                            
                            ident_mat <- diag(lambda, nrow = ncol(x))
                            au <- all.vars(formula)[1]
                            y <- as.matrix(data[,au])
                            beta_ridge_au <- solve(t(R) %*% R + ident_mat) %*% (t(x) %*% y)
                            beta_ridge <<- beta_ridge_au[,1]
                            
                            # yhat
                            
                            y_hat_aux <- x %*% beta_ridge
                            y_hat <<- y_hat_aux[,1]
                            
                          },
                          
                          # print function
                          
                          print = function(){
                            cat("\n","Call:","\n",
                                paste("ridgereg(", "formula = ", formula[2]," ", formula[1], " ", formula[3],", ", "data = ", data_name, ", lambda = ", lambda, ")",sep = "", collapse = "\n" ),
                                "\n","Coefficients:","\n")
                            beta_ridge
                          },
                          
                          # predict function
                          
                          predict = function(df = NULL){
                            
                            results <- y_hat
                            if(!(is.null(df))){
                              df <- data.frame(Intercept = 1, df)
                              a_1 <- as.matrix(df)
                              a_2 <- matrix(beta_ridge, nrow = length(beta_ridge))
                              res <- a_1 %*% a_2
                              result <- res[,1]
                            }
                            return(result)
                          },
                          
                          # Coefficients
                          
                          coef = function(){
                            return(beta_ridge)
                          }
                          
                        )
)
  


