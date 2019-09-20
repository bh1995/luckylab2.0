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
