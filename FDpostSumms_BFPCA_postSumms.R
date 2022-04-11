postSumms <- function(K,                 # Number of eigencomponents, K (integer) 
                      model_components,  # Second element of the list returned from the MCMC 
                                         # function, model_components (list, length 4), containing 
                                         # the M posterior estimates of the model components
                      alpha              # Alpha-level of (1 - alpha)100% credible intervals (scalar, 0-1)
                      ){
  #############################################################################
  ## Description: Function for calculating the traditional and proposed functional 
  ##              depth posterior summaries for the model components described in 
  ##              "Functional Depth Posterior Summaries for Bayesian Functional 
  ##              Principal Component Analysis" by Boland et al. (2022). Calculations 
  ##              included are the traditional posterior summaries for the mean function, 
  ##              covariance function, eigenvalues, and eigenfunctions detailed in 
  ##              Section 2.1 and the proposed posterior summaries for the mean function
  ##              and eigenfunctions detailed in Section 3. Specifically, traditional 
  ##              point estimates are calculated for all the model components, and 
  ##              both traditional and proposed credible intervals are calculated 
  ##              for the mean function and eigenfunctions. For quick reference,
  ##              Table 1 in the manuscript contains a summary of the nomenclature 
  ##              and notation of all the posterior summaries.
  ## Args:        (see above)
  ## Returns:     list()
  ##              mu: (list, length 2)
  ##                point_ests: (list, length 2)
  ##                  mu.hat: estimated mean, \hat{\mu}(t) (vector, T x 1)
  ##                  mu.MBD: MBD median mean, \hat{m}{\mu(t)} (vector, T x 1)
  ##                credible_ints: (list, length 5)
  ##                    pointwise_parametric: P^p_{1 - \alpha}{\mut(t)} (list, length 2)
  ##                      lower: lower bound (vector, T x 1)
  ##                      upper: upper bound (vector, T x 1)
  ##                    simultaneous_parametric: P^s_{1 - \alpha}{\mu(t)} (list, length 2)
  ##                      lower: lower bound (vector, T x 1)
  ##                      upper: upper bound (vector, T x 1)
  ##                    pointwise_quantile: Q^p_{1 - \alpha}{\mut(t)} (list, length 2)
  ##                      lower: lower bound (vector, T x 1)
  ##                      upper: upper bound (vector, T x 1)
  ##                    simultaneous_quantile: Q^s_{1 - \alpha}{\mu(t)} (list, length 2)
  ##                      lower: lower bound (vector, T x 1)
  ##                      upper: upper bound (vector, T x 1)
  ##                    MBD: D_{1 - \alpha}{\mu(t)} (list, length 2)
  ##                      lower: lower bound (vector, T x 1)
  ##                      upper: upper bound (vector, T x 1)
  ##              C: (list, length 1)
  ##                point_ests: (list, length 1)
  ##                  C.hat: estimated covariance, \hat{C}(s, t) (matrix, T x T)
  ##              rho: (list, length K)
  ##                rho_1: (list, length 1)
  ##                ... 
  ##                rho_k: (list, length 1)
  ##                  point_ests: (list, length 2)
  ##                    rho.hat: estimated eigenvalue, \hat{\rho}_k (scalar)
  ##                    rho.tilde: estimated eigenvalue via covariance, \tilde{\rho}_k  (scalar)
  ##                ... 
  ##                rho_K: (list, length 1)
  ##              psi: (list, length K)
  ##                psi_1: (list, length 2)
  ##                ...
  ##                psi_k: (list, length 2)
  ##                  point_ests: (list, length 4)
  ##                    psi.hat: estimated eigenfunction, \hat{\psi}_k(t) (vector, T x 1)
  ##                    psi.tilde: estimated eigenfunction via covariance, \tilde{\psi}_k(t) (vector, T x 1)
  ##                    psi.MBD: MBD median eigenfunction, \hat{m}{\psi_k(t)} (vector, T x 1)
  ##                    psi.MVD: MVD median eigenfunction, \tilde{m}{\psi_k(t)} (vector, T x 1)
  ##                  credible_ints: (list, length 6)
  ##                    pointwise_parametric: P^p_{1 - \alpha}{\psi_k(t)} (list, length 2)
  ##                      lower: lower bound (vector, T x 1)
  ##                      upper: upper bound (vector, T x 1)
  ##                    simultaneous_parametric: P^s_{1 - \alpha}{\psi_k(t)} (list, length 2)
  ##                      lower: lower bound (vector, T x 1)
  ##                      upper: upper bound (vector, T x 1)
  ##                    pointwise_quantile: Q^p_{1 - \alpha}{\psi_k(t)} (list, length 2)
  ##                      lower: lower bound (vector, T x 1)
  ##                      upper: upper bound (vector, T x 1)
  ##                    simultaneous_quantile: Q^s_{1 - \alpha}{\psi_k(t)} (list, length 2)
  ##                      lower: lower bound (vector, T x 1)
  ##                      upper: upper bound (vector, T x 1)
  ##                    MBD: D_{1 - \alpha}{\psi_k(t)} (list, length 2)
  ##                      lower: lower bound (vector, T x 1)
  ##                      upper: upper bound (vector, T x 1)
  ##                    MVD: D^{\star}_{1 - \alpha}{\psi_k(t)} (list, length 2)
  ##                      lower: lower bound (vector, T x 1)
  ##                      upper: upper bound (vector, T x 1)
  ##                ...
  ##                psi_K: (list, length 2)
  ## postSumms Outline:
  ##              1. Traditional posterior summaries
  ##              2. Traditional credible intervals
  ##              3. Proposed point estimates
  ##              4. Proposed credible intervals
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("acid")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.packages) 
  
  # Load packages
  library(acid)
  
  # Create list for output
  postSumms <- list(
    mu = list(point_ests = list(mu.hat = NA, 
                                mu.MBD = NA), 
              credible_ints = list(pointwise_parametric = NA, 
                                   simultaneous_parametric = NA,
                                   pointwise_quantile = NA,
                                   simultaneous_quantile = NA,
                                   MBD = NA)),
    C = list(point_ests = list(C.hat = NA)),
    rho = list(),
    psi = list())
  for(k in 1:K){
    postSumms$rho[[k]] <- list(point_ests = list(rho.hat = NA,
                                                 rho.tilde = NA))
    postSumms$psi[[k]] <- list(
      point_ests = list(psi.hat = NA, 
                        psi.tilde = NA,
                        psi.MBD = NA,
                        psi.MVD = NA), 
      credible_ints = list(pointwise_parametric = NA,
                           simultaneous_parametric = NA,
                           pointwise_quantile = NA,
                           simultaneous_quantile = NA,
                           MBD = NA, MVD = NA))
  }
  names(postSumms$rho) <- paste("rho_", 1:K, sep = "")
  names(postSumms$psi) <- paste("psi_", 1:K, sep = "")
  
  
  #############################################################################
  # 1. Traditional point estimates
  #############################################################################
  
  
  # a. Estimated mean, \hat{\mu}(t) (vector, T x 1)
  mu.hat <- colMeans(model_components$mu.m)
  postSumms$mu$point_ests$mu.hat <- mu.hat
  
  
  # b. Estimated covariance, \hat{C}(s, t) (matrix, T x T)
  C.hat <- Reduce("+", model_components$C.m)/length(model_components$C.m)
  postSumms$C$point_ests$C.hat <- C.hat
  
  
  # c. Eigenvalues
  
  ## i. Estimated eigenvalue, \hat{\rho}_k (scalar)
  for(k in 1:K){
    postSumms$rho[[k]]$point_ests$rho.hat <- mean(model_components$rho.m[, k])
  }
  
  ## ii. Estimated eigenvalue via covariance, \tilde{\rho}_k (scalar)
  SVD.C.hat <- eigen(C.hat)  # Apply SVD to the estimated covariance
  for(k in 1:K){
    postSumms$rho[[k]]$point_ests$rho.tilde <- SVD.C.hat$values[k] * (1/length(mu.hat))
  }
  
  
  # d. Eigenfunctions
  
  ## i. Estimated eigenfunction, \hat{\psi}_k(t) (vector, T x 1)
  psi.hat <- lapply(model_components$psi.m, colMeans)
  for(k in 1:K){
    postSumms$psi[[k]]$point_ests$psi.hat <- psi.hat[[k]]
  }

  ## ii. Estimated eigenfunction via covariance, \tilde{\psi}_k(t) (vector, T x 1)
  for(k in 1:K){
    if(sum((psi.hat[[k]] - SVD.C.hat$vectors[, k])^2)/sum(psi.hat[[k]]^2) >
       sum((psi.hat[[k]] - -1 * SVD.C.hat$vectors[, k])^2)/sum(psi.hat[[k]]^2)){
      SVD.C.hat$vectors[, k] <- -1 * SVD.C.hat$vectors[, k]
    }
    postSumms$psi[[k]]$point_ests$psi.tilde <- SVD.C.hat$vectors[, k]
  }

  #############################################################################
  # 2. Traditional credible intervals
  #############################################################################
  
  # Define function to calculate pointwise parametric credible intervals
  pointwise_parametric <- function(data,  # Sample of posterior estimates, (matrix, M x T)
                                   alpha  # Alpha-level of credible interval, (scalar, 0-1)
                                   ){
    point.mean <- colMeans(data)  # Calculate pointwise mean, (vector, T x 1)
    point.sd <- apply(data, 2, sd)  # Calculate pointwise standard deviation, (vector, T x 1)
    
    # Return list() containing upper (vector, T x 1) and lower (vector, T x 1) bounds
    return(list(lower = point.mean - qnorm(1 - alpha/2) * point.sd,
                upper = point.mean + qnorm(1 - alpha/2) * point.sd))
  }
  
  # Define function to calculate simultaneous parametric credible intervals
  simultaneous_parametric <- function(data,  # Sample of posterior estimates, (matrix, M x T)
                                      alpha  # Alpha-level of credible interval, (scalar, 0-1)
                                      ){
    point.mean <- colMeans(data)  # Calculate pointwise mean, (vector, T x 1)
    point.sd <- apply(data, 2, sd)  # Calculate pointwise standard deviation, (vector, T x 1)
    
    # Calculate sample quantile of the absolute standardized deviation, c_alpha
    c_alpha.m <- unlist(lapply(1:nrow(data), function(m){
      max(abs((data[m,] - point.mean))/point.sd)
    }))
    c_alpha <- quantile(c_alpha.m, probs = 1 - alpha)
    
    # Return list() containing upper (vector, T x 1) and lower (vector, T x 1) bounds
    return(list(lower = point.mean - c_alpha * point.sd,
                upper = point.mean + c_alpha * point.sd))
  }

  
  # a. Mean function
  
  ## i. Parametric credible intervals
  postSumms$mu$credible_ints$pointwise_parametric <- pointwise_parametric( 
    model_components$mu.m, alpha = alpha)
  postSumms$mu$credible_ints$simultaneous_parametric <- simultaneous_parametric(
    model_components$mu.m, alpha = alpha)

  ## ii. Quantile credible intervals
  postSumms$mu$credible_ints$pointwise_quantile <- confband.pw(
    sample = model_components$mu.m, level = 1 - alpha)
  postSumms$mu$credible_ints$simultaneous_quantile <- confband.kneib(
    sample = model_components$mu.m, level = 1 - alpha
  )
  
  
  # b. Eigenfunctions
  
  ## i. Parametric credible intervals
  for(k in 1:K){
    postSumms$psi[[k]]$credible_ints$pointwise_parametric <- pointwise_parametric(
      model_components$psi.m[[k]], alpha = alpha)
    postSumms$psi[[k]]$credible_ints$simultaneous_parametric <- simultaneous_parametric(
      model_components$psi.m[[k]], alpha = alpha)
  }
  
  ## ii. Quantile credible intervals
  for(k in 1:K){
    postSumms$psi[[k]]$credible_ints$pointwise_quantile <- confband.pw(
      sample = model_components$psi.m[[k]], level = 1 - alpha)
    postSumms$psi[[k]]$credible_ints$simultaneous_quantile <- confband.kneib(
      sample = model_components$psi.m[[k]], level = 1 - alpha)
  }
  
  
  #############################################################################
  # 3. Proposed point estimates
  #############################################################################
  ## Note: the following code is adapted from the R code provided in "Exact fast 
  ##       computation of band depth for large functional datasets: How quickly can 
  ##       one million curves be ranked?" Sun et al. (2012)
  
  # Define function to calculate total number of bands that can be formed 
  combinat <- function(n,  # Total number of posterior estimates, M
                       p)  # Total number of curves included in a band (default = 2)
    { # tota
    if (n < p) {combinat = 0}
    else {combinat = exp(lfactorial(n) - (lfactorial(p) + lfactorial(n-p)))}
  }
  
  # Define function to calculate MBD of a sample of curves. 
  MBD <- function(data  # Transposed functional posterior estimates, matrix (T x M)
                  ){
    p = dim(data)[1]  
    n = dim(data)[2]  
    rmat = apply(data, 1, rank)
    down = rmat - 1
    up = n - rmat
    mbd <- (rowSums(up * down) / p + n - 1) / combinat(n, 2)  # Calculate MBD values 
    mbd.rank <- rank(mbd)  # Rank the MBD values from smallest to largest (1 -> M)
    
    # Return data.frame with columns c("id", "mbd", "mbd.rank")
    #        id: row index, 1,...,M
    #        mbd: MBD value, MBD_{M, 2}{X^{(m)}(t)}
    #        mbd.rank: rank of the MBD value, MBD_{M, 2}{X^{[m]}(t)} for function X^{(m)}(t)
    return(data.frame(id = 1:n, mbd = mbd, rank = mbd.rank))
  }
  
  # Define function to calculate MVD of a sample of covariance functions
  MVD <- function(list  # Functional posterior covariance estimates, (list, length M)
                  ){
    data <- do.call(cbind, lapply(1:length(list), function(m){  # Vectorize the covariance surfaces
      as.vector(list[[m]])
    }))
    
    p = dim(data)[1]  
    n = dim(data)[2]  
    rmat = apply(data, 1, rank)
    down = rmat - 1
    up = n - rmat
    mvd <- (rowSums(up * down) / p + n - 1) / combinat(n, 2)  # Calculate MVD values 
    mvd.rank <- rank(mvd)  # Rank the MVD values from smallest to largest (1 -> M)
    
    # Return data.frame with columns c("id", "mvd", "mvd.rank")
    #        id: row index, 1,...,M
    #        mvd: MBD value, MVD_{M, 2}{X^{(m)}(s,t)}
    #        mvd.rank: rank of the MVD value, MVD_{M, 2}{X^{[m]}(s,t)} for function X^{(m)}(s,t)
    return(data.frame(id = 1:n, mvd = mvd, rank = mvd.rank))
  }
  
  
  # a. MBD median mean, \hat{m}{\mu(t)} (vector, T x 1)
  mu.m.MBD <- MBD(t(model_components$mu.m))  # Calculate the MBD
  mu.m.MBD <- mu.m.MBD[order(-mu.m.MBD$rank),]  # Order rows from largest to smallest MBD values 
  postSumms$mu$point_ests$mu.MBD <- model_components$mu.m[mu.m.MBD$id[1],]  # Obtain median
  
  
  # b. Eigenfunctions
  
  ## i. MBD median eigenfunction, \hat{m}{\psi_k(t)} (vector, T x 1)
  psi.m.MBD <- lapply(1:K, function(k){ 
    psi_k.m.MBD <- MBD(t(model_components$psi.m[[k]]))  # Calculate the MBD for each k = 1,...,K 
    psi_k.m.MBD <- psi_k.m.MBD[order(-psi_k.m.MBD$rank),]  # Order rows from largest to smallest MBD values 
    return(psi_k.m.MBD)
  })
  for(k in 1:K){
    postSumms$psi[[k]]$point_ests$psi.MBD <- model_components$psi.m[[k]][psi.m.MBD[[k]]$id[1],] # Obtain median
  }
  
  ## ii. MVD median eigenfunction, \tilde{m}{\psi_k(t)} (vector, T x 1)
  C.m.MVD <- MVD(model_components$C.m)  # Calculate the MVD of the covariance functions
  C.m.MVD <- C.m.MVD[order(-C.m.MVD$rank),]  # Order rows from largest to smallest MVD values 
  for(k in 1:K){
    postSumms$psi[[k]]$point_ests$psi.MVD <- model_components$psi.m[[k]][C.m.MVD$id[1],] # Obtain median
  }
  
  #############################################################################
  # 4. Proposed credible intervals
  #############################################################################
  
  # Define function to calculate MVD or MBD credible interval
  depth_credible_interval <- function(rank,  # Data.frame returned from the MBD() or MVD() functions
                                      data,  # Posterior estimates (matrix, M x T)
                                      alpha  # Alpha-level of credible intervals
                                      ){
    
    M <- nrow(data)  # Number of posterior estimates, M
    
    # Obtain the row indices for the posterior estimates with the largest MBD/MVD values
    row.indices <- rank$id[rank$rank %in% ((alpha * M) + 1):M]
    
    # Subset of posterior estimates with largest MBD or MVD to form credible interval 
    data.subset <- data[row.indices, ]
    
    # Return list of lower (T x 1) and upper (T x 1) bounds
    return(list(lower = apply(data.subset, 2, min),
                upper = apply(data.subset, 2, max)))
  }
  
  
  # a. Mean function MBD credible interval
  postSumms$mu$credible_ints$MBD <- depth_credible_interval(
    mu.m.MBD, model_components$mu.m, alpha)
  
  
  # b. Eigenfunctions
  for(k in 1:K){
    
    ## i. Calculate MBD credible interval
    postSumms$psi[[k]]$credible_ints$MBD <- depth_credible_interval(
      rank = psi.m.MBD[[k]], model_components$psi.m[[k]], alpha)
    
    ## ii. Calculate MVD credible interval
    postSumms$psi[[k]]$credible_ints$MVD <- depth_credible_interval(
      C.m.MVD, model_components$psi.m[[k]], alpha)
  }
  
  return(postSumms)
  
}