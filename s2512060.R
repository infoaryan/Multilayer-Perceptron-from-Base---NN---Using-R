# Name: Aryan Verma, Student  Number: s2512060

# Function to compute η(r)
# Input - number
# Output - number
eta <- function(r) {
  if (r > 0) {
    return(r^2 * log(r))
  } else {
    return(0)
  }
}


# getTPS function
# Input: x (n × 2 matrix) and k (No. of basis functions)
# Output: xk (k × 2 matrix), X matrix, S matrix, Z Matrix
getTPS <- function(x, k = 100) {
  
  # getting number of rows in x
  n <- nrow(x)
  
  # Resetting k to n if it is greater than or equal to n
  if (k >= n) {
    k <- n
    x_star <- x
  } else {
    # Select k rows randomly from x
    sampled_rows <- sample(1:n, k)
    x_star <- x[sampled_rows, ]
  }
  
  # Compute T-(1,x) and Ts-(1,x_star) matrix
  Ts <- cbind(1, x_star)
  T_ <- cbind(1, x)
  
  # Compute Z - the last k−3 columns of Q from QR Decomposition of Ts
  Z <- qr.Q(qr(Ts),complete=TRUE)[,-(1:3)]
  
  # Calculate E matrix, E=η(x), where η(x) = η(||x − x∗j||).
  E <- matrix(0, nrow(x), nrow(x_star))
  for (i in 1:nrow(x)) {
    for (j in 1:nrow(x_star)) {
      r <- sqrt((x[i, 1] - x_star[j, 1])^2 + (x[i, 2] - x_star[j, 2])^2)
      E[i, j] <- eta(r)
    }
  }
  
  # Calculate E* Matrix, E=η(x), where η(x) = η(||x − x∗j||)
  E_star <- matrix(0, nrow(x_star), nrow(x_star))
  for (i in 1:nrow(x_star)) {
    for (j in 1:nrow(x_star)) {
      r <- sqrt((x_star[i, 1] - x_star[j, 1])^2 + (x_star[i, 2] - x_star[j, 2])^2)
      E_star[i, j] <- eta(r)
    }
  }
  
  # Combining E*Z and T to create the X matrix
  X <- cbind(E %*% Z, T_)
  
  # Computing QR decomposition of X
  qr_result <- qr(X)
  
  # Extract Q from QR decomposition
  Q <- qr.Q(qr_result)
  
  # Compute S matrix
  S <- t(Z) %*% E_star %*% Z
  # Padding the dimensions for zeros
  res_row = k-nrow(S)
  res_col = k-ncol(S)
  S <- cbind(S, matrix(0,nrow(S), res_col))
  S<- rbind(S,matrix(0, res_row, ncol(S)))
  
  # Returning the named list of required elements
  result <- list(
    xk = x_star,
    X = X,
    S = S,
    Z = Z
  )
  return(result)
}


# fitTPS function
# INPUT: x (n × 2 matrix of location values)
#        y (n vector of values to smooth)
#        k - No. of basis functions to use
#        lsp - log λ limits for smoothening param
# OUTPUT: An object of class tps, with list of named items:
#         beta - best fit β at the optimal λ value
#         mu - corresponding μ
#         medf - effective degrees of freedom at the optimal λ 
#         lambda - vector of 100 λ values
#         gcv - corresponding GCV score vector
#         edf - corresponding vector of EDFs
fitTPS <- function(x,y,k=100,lsp=c(-5,5)){

  # initialising the required variables and vectors
  min_gcv <- Inf
  best_lambda <- NULL
  best_beta <- NULL
  best_mu <- NULL
  best_edf <- NULL
  gcv <- numeric(100)
  TPS <- getTPS(x,k)
  EDF_vec <- numeric(100)
  
  # Getting the possible values of lambda in log range of lsp
  lambda_grid <- 10^seq(lsp[1], lsp[2], length = 100)
  
  # QR Decompose X 
  qr_decom = qr(TPS$X)
  
  # Extract Q and R matrix from the QR decomposition
  R <- qr.R(qr_decom)
  Q <- qr.Q(qr_decom)
  
  # Perform the symmetric eigen decomposition
  evd <- eigen(t(solve(R)) %*% TPS$S %*% solve(R))
  
  # Extract U, Λ, and U^T from the eigen decomposition
  U <- evd$vectors
  Λ <- diag(evd$values)
  Ut <- t(U)
  
  # Finding optimum lambda value by looping through all possible values
  for (i in 1:length(lambda_grid)) {

    lambda = lambda_grid[i]
    
    # Calculating the EDF (Effective degrees of freedom)
    EDF <- t(sum(1 / (1 + lambda * diag(Λ))))
    
    # Calculating beta value
    # Calculate the matrix (I + λΛ)^(−1)
    temp1 <- solve(diag(1, ncol(Λ)) + lambda * Λ)
    # Calculate the matrix U * (I + λΛ)^(−1) * U^T
    temp2 <- U %*% temp1 %*% Ut
    # Calculate Q^T * y
    Qt_y <- t(Q) %*% y
    # Calculate R^(-1) * (U * (I + λΛ)^(−1) * U^T) * Q^T * y
    beta_hat <- solve(R) %*% (temp2 %*% Qt_y)
    
    # Calculating the mu - model prediction
    mu <- TPS$X %*% beta_hat
    
    #Calculating the gcv_score
    gcv_score  <- sum((y - mu)^2) / (length(y) - EDF)^2
    gcv_score <- Re(gcv_score)
    gcv[i] <- gcv_score
    EDF_vec[i] <- EDF
    
    
    # Update the minimum GCV and best_lambda if a lower GCV value is found
    if (gcv_score < min_gcv) {
      min_gcv <- gcv_score
      best_lambda <- lambda
    }
  }
  
  ############################
  ### Calulating best_beta ###
  ############################
  temp1 <- solve(diag(1, ncol(Λ)) + best_lambda * Λ)
  # Calculate the matrix U * (I + λΛ)^(−1) * U^T
  temp2 <- U %*% temp1 %*% Ut
  # Calculate Q^T * y
  Qt_y <- t(Q) %*% y
  # Calculate R^(-1) * (U * (I + λΛ)^(−1) * U^T) * Q^T * y
  best_beta <- solve(R) %*% (temp2 %*% Qt_y)
  
  ############################
  #### Calulating best_mu ####
  ############################
  best_mu <- TPS$X %*% best_beta
  
  
  ############################
  ### Calulating best_medf ###
  ############################
  best_edf <- t(sum(1 / (1 + best_lambda * diag(Λ))))
  
  tps_obj <- list(
    beta = best_beta,
    mu = best_mu,
    medf = best_edf,
    lambda = lambda_grid,
    gcv = gcv,
    edf = EDF_vec,
    Z = TPS$Z,
    xk = TPS$xk
  )
  class(tps_obj) <- "tps"  # Set the class of the object to "tps"
  return(tps_obj)
}
  

# plot.tps function
# INPUT: tps_obj - An object of class tps
#        m - No. of points 
#.       lower_sequence - lower limit for generation of new sequential data
#.       upper_sequence - lower limit for generation of new sequential data
# OUTPUT: Perspective plot of Fitted Thin Plate Spline
plot.tps <- function(tps_obj, m=50, lower_sequence=0, upper_sequence=1) {
  
  # Extract values from the TPS object
  beta <- tps_obj$beta
  xk <- tps_obj$xk
  
  # Now, extract the β(beta) to α(alpha) and δ(delta)
  cutoff = nrow(beta)
  # Select last 3 values to be alpha
  alpha <-beta[(cutoff-2):cutoff,]
  # Rest of the values are put to detla_cap
  delta_cap <- beta[1:(cutoff-3)]
  # Now, delta is calculated from delta_cap
  delta <- tps_obj$Z %*% delta_cap
  
  # Defining the new data to be visualized
  x2 <- x1 <- seq(lower_sequence,upper_sequence,length=m)
  xp <- cbind(rep(x1,m),rep(x2,each=m))
  
  # Calculating the predictions on new data
  T_new <- cbind(1,xp)
  
  #Calculate E matrix
  E_new <- matrix(0, nrow(xp), nrow(xk))
  for (i in 1:nrow(xp)) {
    for (j in 1:nrow(xk)) {
      r <- sqrt((xp[i, 1] - xk[j, 1])^2 + (xp[i, 2] - xk[j, 2])^2)
      E_new[i, j] <- eta(r)
    }
  }
  # Prediction T*alpha + E*delta
  prediction <- T_new %*% alpha+ E_new %*% delta
  # Taking the Real part of prediction only 
  # There are some times values with 0 imaginary part
  # We are just taking real part from those predicted values
  prediction <- Re(prediction)
  
  # Plotting the perspective plot 
  #contour(x1,x2,matrix(prediction, m, m))
  persp(x1,x2,matrix(prediction,m,m),theta=30,phi=30)
  
}


# INPUT: None
# OUTPUT: perspective plot on TPS using sample data 
# This function will use the sample data to run through the code
testExample <- function(){
  # Generatign sample data for fitting the spline
  ff <- function(x) exp(-(x[,1]-.3)^2/.2^2-(x[,2] - .3)^2/.3^2)*.5 + exp(-(x[,1]-.7)^2/.25^2 - (x[,2] - .8 )^2/.3^2)
  n <- 500
  x <- matrix(runif(n*2),n,2)
  y <- ff(x) + rnorm(n)*.1 
  
  # Fitting the TPS
  tps_result <- fitTPS(x,y)
  
  # Plotting the fitted thin plate spline
  plot.tps(tps_result)
}

testExample()

