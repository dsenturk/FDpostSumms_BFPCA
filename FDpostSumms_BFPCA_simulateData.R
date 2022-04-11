simulateData <- function(n, # number of subjects (scalar)
                         t  # functional domain grid (T x 1 vector)
                         ){
  #############################################################################
  ## Description: Function for generating densely observed functional data as described 
  ##              in Supporting Information Appendix D. In particular, the simulated 
  ##              data contain no outliers (Case 1) and are generated from the standard 
  ##              FPCA model detailed in Section 2.1.
  ## Args:        (see above)
  ## Returns:     Y: matrix of simulated functional data, Y_i(t_j) (matrix, n x T)
  ## simulateData Outline:
  ##              1. Define model components
  ##              2. Generate functional data
  #############################################################################

  #############################################################################
  # 1. Define model components
  #############################################################################
  
  # a. Global variables
  K <- 2  # Number of eigencomponents 
  T_tot <- length(t)  # Number of time points, T
  
  # b. Mean function, mu(t) (vector, T x 1)
  mu_vec <- 10 * sqrt(1 - 2 * (t - 0.5)^2)
  
  # c. Eigenfunctions, psi_k(t) (vector, T x 1), and eigenvalues, rho_k (scalar)
  psi_1 <- sqrt(2) * sin(2 * pi * t)
  rho_1 <- 15 
  psi_2 <- sqrt(2) * cos(2 * pi * t)
  rho_2 <- 5  
  psi <- cbind(psi_1, psi_2)  # Matrix of eigenfunctions (matrix, T x K)
  rho <- diag(c(rho_1, rho_2))  # Diagonal matrix of eigenvalues (matrix, K x K)
  
  # d. Measurement error, sigma^2_{epsilon} (scalar)
  sig2_e <- 15
  
  #############################################################################
  # 2. Generate functional data
  #############################################################################
  
  # a. Simulate FPCA scores, \xi_{ik} (matrix, n x K)
  xi_i <- matrix(rnorm(n = n * K, mean = 0, sd = 1), ncol = K) %*% sqrt(rho)
  
  # b. Calculate smooth underlying functions, f_i(t) (matrix, n x T)
  f_i <- xi_i %*% t(psi) + t(matrix(rep(mu_vec, n), nrow = T_tot)) 
  
  # c. Generate random noise and construct functional observations, Y_i(t) (matrix, n x T)
  Y <- do.call(rbind, lapply(1:n, function(x) f_i[x, ] + sqrt(sig2_e) * rnorm(T_tot)))
  
  return(Y)
}