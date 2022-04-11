MCMC <- function(Y,           # Matrix of functional observations Y_i(t_j) (matrix, n x T)
                 t,           # Functional domain grid t_j, j = 1,...,T (vector, T x 1)
                 runs,        # Total number of iterations in MCMC chain (integer)
                 burnin,      # Percentage of total iterations used for burn-in (scalar, 0-1)
                 thin,        # Thinning interval between consecutive iterations (integer)
                 r_tot = NA,  # Number of b-spline basis functions, R (integer), automatically 
                              # chosen as R = floor(T/2) if value is left as NA
                 l_tot = NA   # Number of latent components, L (integer), automatically 
                              # chosen as L = max{floor(R/2), 6} if value is left as NA
){
  #############################################################################
  ## Description: Function for performing posterior estimation using a Gibbs sampler 
  ##              described in "Functional Depth Posterior Summaries for Bayesian 
  ##              Functional Principal Component Analysis" by Boland et al. (2022). 
  ##              The sampling model, prior distributions, and posterior distributions
  ##              employed and calculated are described in Section 2.1 and Appendix 
  ##              A of the manuscript. Posterior samples of size M using the BFPCA 
  ##              model are returned from this function and include the set of model 
  ##              parameters and model components described in section 2.1 as well 
  ##              as the post-processing posterior alignment of the eigenfunctions 
  ##              detailed in Appendix B. 
  ## Args:        see above
  ## Returns:     list()
  ##              model_params: (list, length 9)
  ##                beta.m: mean vector coefficients, \beta^{(m)} (matrix, M x R) 
  ##                Lambda.m: vectorized factor loading values, Vec(\Lambda^{(m)}) (matrix, M x RL)
  ##                eta.m: vectorized subject-specific scores, Vec(\eta_1^{(m)},...,\eta_n^{(m)}) (matrix, M x nL)
  ##                sig2_e.m: measurement error variance, \sigma^{2(m)}_{\epsilon} (vector, M x 1)
  ##                sig2_beta.m: variance component for mean vector coefficients, \sigma^{2(m)}_{\beta} (vector, M x 1)
  ##                sig2_lambda_rl.m: variance components for factor loadings, \sigma^{2(m)}_{\varphi_{r \ell}} (matrix, M x RL)
  ##                varphi.m: factor loading element-wise precisions, \varphi_{r \ell} (matrix, M x RL)
  ##                tau.m: factor loading column-wise precision, \tau_{\ell} (matrix, M x L)
  ##                delta.m: multiplicative gamma process, \delta_h (matrix, M x L)
  ##              model_components: (list, length 4)
  ##                mu.m: mean function estimates, \mu^{(m)}(t) (matrix, M x T) 
  ##                C.m: covariance function estimates, C^{(m)}(s,t) (list, length M) 
  ##                  C^{(1)}(s,t): (matrix, T x T)
  ##                  ...
  ##                  C^{(M)}(s,t): (matrix, T x T)
  ##                rho.m: eigenvalue estimates, \rho^{(m)} (matrix, M x K) 
  ##                psi.m: eigenfunction estimates, \psi^{(m)}_k(t) (list, length K) 
  ##                  psi_1: estimates for \psi^{(m)}_1(t) (matrix, M x T)
  ##                  ...
  ##                  psi_K: estimates for \psi^{(m)}_K(t) (matrix, M x T)
  ## MCMC Outline:
  ##              1. Define global variables and hyperparameters
  ##              2. Initialize MCMC Chain values 
  ##              3. Employ Gibbs sampler 
  ##              4. Estimate posterior model components
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("fda", "MASS")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.packages) 
  
  # Load packages
  library(fda)
  library(MASS)
  
  # Define function for drawing from truncated gamma distribution
  rgammat <- function(n, range, shape, scale = 1) {
    
    F.a <- pgamma(min(range), shape = shape, scale = scale)
    F.b <- pgamma(max(range), shape = shape, scale = scale)
    
    u <- runif(n, min = F.a, max = F.b)
    qgamma(u, shape = shape, scale = scale)
  }
  
  #############################################################################
  # 1. Define global variables and hyperparameters
  #############################################################################

  # a. Define global variables
  i_tot <- nrow(Y)  # Total number of subjects, n
  j_tot <- ncol(Y)  # Total number of time points, T
  if(is.na(r_tot)){
    r_tot <- floor(j_tot/2) # Total number of basis functions, R
  }
  if(is.na(l_tot)){
    l_tot <- max(6, floor(r_tot/4))  # Total number of latent components, L
  }
  
  
  # b. Define posterior sampling index after burn-in and thinning
  burnin_iters <- burnin * runs  # Last burn-in iteration (scalar)
  m <- seq(from = burnin_iters + 1, to = runs, by = thin)  # Posterior sample index (M x 1)
  m_tot <- length(m)  # Total number of posterior estimates after burn-in and thinning (scalar)
  
  
  # c. Create matrix of B-spline basis functions, B (T x R)
  B <- unname(eval.basis(t, create.bspline.basis(rangeval = c(min(t), max(t)),
                                                 nbasis = r_tot))) 
  
  
  # d. Define prior hyperparameters
  
  ## i. Mean coefficient vector prior, \beta (R x 1)
  Omega <- bandSparse(n = r_tot, k = -1:0, symm = TRUE,  # positive definite (R x R) penalty matrix
                      diag = list(rep(-1, r_tot - 1), 
                                  c(1, rep(2, r_tot - 2), 1))) + 
    Diagonal(n = r_tot, x = 0.00001)
  Omega <- as.matrix(Omega) 
  a_beta <- 2.0  # Variance component hyperparameter, \sigma^2_{\beta}
  
  ## ii. Loading factor element-wise precision, \varphi_{r \ell} 
  nu <- 10  # Variance component hyperparameter
  
  ## iii. Loading factor column-wide precision, tau_{\ell}, multiplicative gamma distributions
  a_1 <- 1.1  # Variance component hyperparameter, \delta_1
  a_2 <- 2  # Variance component hyperparameter, \delta_h, h >= 2
  
  #############################################################################
  # 2. Initialize MCMC Chain values
  #############################################################################
  
  # a. Variance component, \sigma^2_{\beta} (scalar)
  sig2_beta <- 1/i_tot
  
  
  # b. Mean coefficient vector, \beta (R x 1)
  BtB <- crossprod(B) 
  Ybar <- colMeans(Y)
  BtYbar <- crossprod(B, Ybar)
  beta <- solve(BtB + ((sig2_beta / i_tot) * Omega)) %*% BtYbar
  
  
  # c. Loading factor matrix, \Lambda (R x L)
  Bbeta <- B %*% beta
  Yc <- Y - tcrossprod(rep(1, i_tot), Bbeta)
  theta <- Yc %*% B %*% solve(BtB + (sig2_beta * Omega))
  V_theta <- crossprod(theta)/i_tot
  SVD.V_theta <- eigen(V_theta)
  Lambda <- SVD.V_theta$vectors[, 1:l_tot] %*% 
    sqrt(diag(SVD.V_theta$values)[1:l_tot, 1:l_tot])
  
  
  # d. Measurement error Variance, \sigma^2_{\epsilon} (scalar)
  sig2_e <- sum(diag(crossprod(Yc - tcrossprod(theta, B))))/(i_tot * j_tot)
  
  
  # e. Subject-specific scores, \eta_i (N x L)
  eta <- matrix(rnorm(n = i_tot * l_tot, mean = 0, sd = 1), nrow = i_tot, 
                ncol = l_tot)
  
  
  # f. Factor loading precisions, \varphi_{r \ell} (R x L)
  varphi <- matrix(rgamma(n = r_tot * l_tot, shape = nu/2, scale = 2/nu), 
                   nrow = r_tot, ncol = l_tot)
  
  
  # g. Column-wise precisions, \tau_{\ell}, (L x 1)
  delta <- c(rgamma(n = 1, shape = a_1, scale = 1),
             rgammat(n = l_tot - 1, range = c(1, Inf), shape = a_2, scale = 1))
  tau <- cumprod(delta)
  
  
  # h. Initialize matrices and vectors for storing the M posterior estimates 
  beta.m <- matrix(NA, nrow = m_tot, ncol = r_tot)  # (M x R)
  Lambda.m <- matrix(NA, nrow = m_tot, ncol = r_tot * l_tot)  # (M x RL) vectorized \Lambda per row
  sig2_e.m <- numeric(m_tot)  # (M x 1)
  eta.m <- matrix(NA, nrow = m_tot, ncol = i_tot * l_tot)  # (M x nL) vectorized \eta_i per row
  sig2_beta.m <- numeric(m_tot)  # (M x 1)
  varphi.m <- matrix(NA, nrow = m_tot, ncol = r_tot * l_tot)  # (M x RL) vectorized \varphi_{r \ell} per row
  delta.m <- matrix(NA, nrow = m_tot, ncol = l_tot)  # (M x L)
  
  #############################################################################
  # 3. Employ Gibbs sampler (Appendix B)
  #############################################################################

  for (run in 1:runs) {
    
    #### Step 1. Update mean coefficients, \beta (R x 1) ####
    BLambda <- B %*% Lambda
    Sigma_y <- tcrossprod(BLambda) + (sig2_e * diag(j_tot))
    Sigma_y.inv <- solve(Sigma_y)
    BtSigma_y.inv <- crossprod(B, Sigma_y.inv)
    V_beta.post <- solve(BtSigma_y.inv %*% B +  # Posterior Variance v^{post}_{\beta}
                           (Omega/(i_tot * sig2_beta)))/i_tot  
    mu_beta.post <- (i_tot * V_beta.post) %*% BtSigma_y.inv %*% Ybar  # Posterior mean \mu^{post}_{\beta}
    beta <- mvrnorm(n = 1, mu = mu_beta.post, Sigma = V_beta.post)  # Update value 
    
    
    #### Step 2. Update loading matrix, \Lambda (R x L) ####
    Bbeta <- B %*% beta
    Yc <- Y - tcrossprod(rep(1, i_tot), Bbeta)
    Tau <- matrix(tcrossprod(rep(1, r_tot), tau), nrow = r_tot, ncol = l_tot)
    Sigma_Lambda <- diag(1/as.vector(varphi * Tau)) + 
      0.0000001 * diag(r_tot * l_tot)
    Sigma_Lambda.inv <- solve(Sigma_Lambda)
    etateta <- crossprod(eta)
    V_Lambda.post <- sig2_e * solve(kronecker(etateta, BtB) +  # Posterior Variance v^{post}_{\Lambda}  
                                      sig2_e * Sigma_Lambda.inv)
    mu_Lambda.post <- tcrossprod(V_Lambda.post / sig2_e,  # Posterior mean \mu^{post}_{\Lambda}   
                                 crossprod(as.vector(t(Yc)), kronecker(eta, B)))
    Lambda <- matrix(mvrnorm(n = 1, mu = mu_Lambda.post, Sigma = V_Lambda.post),  # Update value  
                     nrow = r_tot, ncol = l_tot)
    
    
    #### Step 3. Update measurement error variance, \sigma^2_{\epsilon} ####
    BLambda <- B %*% Lambda
    RSS <- sum(diag(crossprod(Yc - tcrossprod(eta, BLambda))))
    a_sig2_e.post <- (i_tot * j_tot) / 2  # Posterior shape parameter, a^{post}_{\sigma^2_{\epsilon}} 
    b_sig2_e.post <- RSS / 2  # Posterior scale parameter, b^{post}_{\sigma^2_{\epsilon}}
    H_e <- rgamma(n = 1, shape = a_sig2_e.post, scale = 1 / b_sig2_e.post)
    sig2_e <- 1 / H_e  # Update value
    
    
    #### Step 4. Update subject-specific scores, \eta_i (N x L) ####
    V_eta.post <- sig2_e * solve(crossprod(BLambda) + sig2_e * diag(l_tot))  # Posterior Variance v^{post}_{\eta_i}   
    mu_eta.post <- Yc %*% BLambda %*% (V_eta.post / sig2_e)  # Posterior mean \mu^{post}_{\eta_i}   
    eta <- t(vapply(X = 1:i_tot, FUN = function(X){  # Update value
      mvrnorm(n = 1, mu = mu_eta.post[X, ], Sigma = V_eta.post)}, 
      numeric(l_tot)))
    
    
    #### Step 5. Update variance component for beta, \sigma^2_{\beta} ####
    a_sig2_beta.post <- (r_tot + a_beta) / 2  # Posterior shape parameter, a^{post}_{\sigma^2_{\beta}}
    b_sig2_beta.post <- (a_beta + crossprod(beta, Omega) %*% beta) / 2  # Posterior scale parameter, b^{post}_{\sigma^2_{\beta}}
    sig2_beta <- rgamma(n = 1, shape = a_sig2_beta.post,  # Update value 
                        scale = 1 / b_sig2_beta.post)
    
    
    #### Step 6. Update factor loading element-wise precisions, varphi_{r \ell} (R x L) ####
    Lambda2 <- Lambda ^ 2
    a_varphi_rl.post <- (nu + 1) / 2  # Posterior shape parameter, a^{post}_{varphi_{r \ell}}
    b_varphi_rl.post <- (nu + (Lambda2 * Tau)) / 2  # Posterior scale parameter, b^{post}_{varphi_{r \ell}}
    varphi <- matrix(rgamma(n = r_tot * l_tot, shape = a_varphi_rl.post,  # Update value 
                            scale = 1 / as.vector(b_varphi_rl.post)),
                  nrow = r_tot, ncol = l_tot, byrow = FALSE)
    
    
    #### Step 7. Update factor loading column-wise precisions, tau_{ell} ####
    
    #### (a) \delta_1
    gamma <- colSums(Lambda2 * varphi)  # \gamma_{\ell}
    a_delta_1.post <- a_1 + ((r_tot * l_tot) / 2)  # Posterior shape parameter, a^{post}_{\delta_1}
    b_delta_1.post <- 1 + (0.5 * (1/delta[1]) * sum(t(tau) * gamma))  # Posterior scale parameter, b^{post}_{\delta_1}
    delta[1] <- rgamma(n = 1, shape = a_delta_1.post,  # Update value 
                       scale = 1 / b_delta_1.post)
    tau <- cumprod(delta)  # Update value 
    
    #### (b) \delta_h, h >= 2
    a_delta_h.post <- a_2 + ((r_tot * (l_tot - 2:l_tot + 1)) / 2) # Posterior shape parameter, a^{post}_{\delta_h}
    for(h in 2:l_tot){
      temp <- t(tau) * gamma
      b_delta_h.post <- 1 + .5 * (1/delta[h]) * sum(temp[h:l_tot])  # Posterior scale parameter, b^{post}_{\delta_h}
      delta[h] <- rgammat(n = 1, range = c(1, Inf),  # Update value  
                          shape = a_delta_h.post[h - 1], 
                          scale = 1 / b_delta_h.post)
      tau <- cumprod(delta)  # Update value 
    }

    
    # Save posterior estimates included after burn-in and thinning
    if(any(run == m)){
      row <- ((run - burnin_iters - 1)/thin) + 1
      beta.m[row, ] <- beta
      Lambda.m[row, ] <- as.vector(Lambda)
      sig2_e.m[row] <- sig2_e
      eta.m[row, ] <- as.vector(eta)
      sig2_beta.m[row] <- sig2_beta
      varphi.m[row,] <- as.vector(varphi)
      delta.m[row,] <- delta
    }
  }
  
  
  # Create list of model parameters for output
  tau.m <- t(apply(delta.m, 1, cumprod))  # \tau_{\ell} (matrix, M x L)
  sig2_lambda_rl.m <- 1/(varphi.m * t(apply(tau.m, 1, rep, each = r_tot)))  # \sigma^2_{\lambda_{r \ell}}(matrix, M x RL)
  model_params <- list(beta.m = beta.m, 
                       Lambda.m = Lambda.m, 
                       eta.m = eta.m,
                       sig2_e.m = sig2_e.m, 
                       sig2_beta.m = sig2_beta.m,
                       sig2_lambda_rl.m = sig2_lambda_rl.m, 
                       varphi.m = varphi.m,
                       tau.m = tau.m, 
                       delta.m = delta.m)
  
  
  print("Gibbs sampling completed")
  
  #############################################################################
  # 4. Estimate posterior model components
  #############################################################################

  # a. Calculate mean function estimates, \mu^{(m)}(t) (matrix, M x T)
  mu.m <- beta.m %*% t(B)  # (matrix, M x T)
  
  
  # b. Calculate covariance function estimates C^{(m)}(s,t), (list, length M)
  C.m <- lapply(1:m_tot, function(y){
    Lambda.i <- matrix(Lambda.m[y,], nrow = r_tot, ncol = l_tot)
    Covfi <- tcrossprod(B %*% Lambda.i)  # (matrix, T x T)
    return(Covfi)  
  })
  names(C.m) <- paste("C^{(", 1:m_tot, ")}(s,t)", sep = "")
  
  
  # c. Decompose covariance function using SVD and obtain eigencomponents 
  
  ## Note: By default obtain K = L eigencomponents as the rank of C(s,t) is <= L,
  ##       but can be modified as needed by the user by replacing l_tot with another 
  ##       value in following code:
  SVD.C.m <- lapply(1:m_tot, function(m){  # (list, length 3)
    C.eigen <- eigen(C.m[[m]])
    rho <- C.eigen$values  
    FVE <- rho/sum(rho)  # fraction of variance explained
    cumFVE <- cumsum(FVE)  
    K <- min(which(cumFVE >= 0.80))  # Number of eigencomponents that explain 80% of the variation
    rho <- rho[1:l_tot] * (1/j_tot)  # (vector, L x 1)
    psi <- C.eigen$vectors[, 1:l_tot]  # (matrix, T x L)
    return(list(K = K, rho = rho, psi = psi))
  })
  
  ## Note: Obtain number of eigencomponents explaining at least 80% of variation in the data 
  ##       for each iteration of the Gibbs sampler if needed for empirical selection by user
  K.m <- unlist(lapply(1:m_tot, function(m){
    SVD.C.m[[m]]$K
  }))
  K.mean <- mean(K.m); K.median <- median(K.m)
  
  ## Note: User can modify K based on empirical considerations, but is chosen to be the median K
  ##       across the M iterations explaining 80% of the variation in our implementation 
  K <- K.median
  
  
  # d. Obtain set of eigenvalues, rho_k (matrix, M x K) 
  rho.m <- do.call(rbind, lapply(1:m_tot, function(m){
    SVD.C.m[[m]]$rho[1:K]}))
  

  # e. Align posterior estimates and obtain set of eigenfunctions, psi_k(t) (list, length K)
  
  ## Define function for posterior alignment of the eigenfunctions (Appendix B)
  eigenAlign <- function(psi_k){  # psi_k: posterior estimates of kth eigenfunction (matrix, M x T) 
    a <- numeric(m_tot)
    a[1] <- 1  # a^{(1)}
    psi_k.bar <- psi_k[1, ] # ergodic mean, \bar{\psi}^{(1)}_k(t)
    d_plus <- sum(abs(psi_k.bar - psi_k[2, ])) # d^{(2)}_{+}
    d_minus <- sum(abs(psi_k.bar + psi_k[2, ])) # d^{(2)}_{-}
    a[2] <- (d_plus <= d_minus) - (d_plus > d_minus)  # a^{(2)}
    psi_k[2,] <- a[2] * psi_k[2,] # aligned estimate, \psi^{(2)\star}_k(t)
    for(m in 3:m_tot){
      psi_k.bar <- colMeans(psi_k[1:(m - 1),])  # ergodic mean, \bar{\psi}^{(m - 1)}_k(t)
      d_plus <- sum(abs(psi_k.bar - psi_k[m, ]))  # d^{(m)}_{+}
      d_minus <- sum(abs(psi_k.bar + psi_k[m, ]))  # d^{(m)}_{-}
      a[m] <- (d_plus <= d_minus) - (d_plus > d_minus)  # a^{(m)}
      psi_k[m,] <- a[m] * psi_k[m,]  # aligned estimate, \psi^{(m)\star}_k(t)
    }
    return(psi_k)  # return: aligned posterior estimates of the kth eigenfunction (matrix, M x T)
  }
  
  ## Obtain posterior eigenfunction samples used for posterior summaries
  psi.m <- lapply(1:K, function(k){
    psi_k <- do.call(rbind, lapply(1:m_tot, function(m){  
      SVD.C.m[[m]]$psi[, k]}))
    psi_k <- eigenAlign(psi_k) # estimates for k (matrix, M x T)
    return(psi_k)
  })
  names(psi.m) <- paste("psi", 1:K, sep = "_")
  
  
  # Create list of model components for output
  model_components <- list(mu.m = mu.m, 
                           C.m = C.m, 
                           rho.m = rho.m, 
                           psi.m = psi.m)
  
  return(list(model_params = model_params, model_components = model_components))
}