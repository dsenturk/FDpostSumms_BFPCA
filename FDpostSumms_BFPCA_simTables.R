simTables <- function(postEsts,     # list() containing the M posterior estimates of the model components across the
                                    # simulation runs where the total number of elements is the number of simulation
                                    # runs (10 in the Tutorial), and each element is the list named model_components 
                                    # (list, length 4) returned from the MCMC() function for the given simulation run.
                      postSumms,    # list() containing the posterior summaries of the model components across the
                                    # simulation runs where each element is the the list returned from the postSumms()
                                    # function for the given run, and the total number of elements is the number of 
                                    # simulation runs (10 in the Tutorial)
                      tableIndex,   # Table index (1 - ISME/MSE, 2 - AR/NAR)
                      truth_funcs,  # list() containing "true" values of model components of the standard FPCA model 
                                    # used in the simData() function when generating the functional data (list, length 3)
                                        # mu: mean function (vector, T x 1)
                                        # rho: (list, length K)
                                          # rho_1 (scalar)
                                          # ...
                                          # rho_K (scalar)
                                        # psi: (list, length K)
                                          # psi_1 (vector, T x 1)
                                          # ...
                                          # psi_K (vector, T x 1)
                      K             # Total number of eigencomponents (integer)
                      ){
  #############################################################################
  ## Description: Function for generating tables containing the ISME/MSE (Table 2) 
  ##              and AR/NAR (Table 3) performance measures to evaluate the performance 
  ##              of the point estimates and credible intervals, respectively. The tables 
  ##              returned from this function correspond to simulation case 1 in Section 
  ##              4 of "Functional Depth Posterior Summaries for Bayesian Functional 
  ##              Principal Component Analysis" by Boland et al. (2022). Function returns
  ##              either the ISME/MSE or the AR/NAR based on the tableIndex argument. 
  ##              table. 
  ## Args:        (see above)
  ## Returns:     list()
  ## simTables Outline:
  ##              1. ISME and MSE Table for Point Estimates
  ##              2. AR and NAR Table for Credible Intervals
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("kableExtra", "tidyverse")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.packages) 
  
  # Load packages
  library(kableExtra)
  library(tidyverse)
  
  #############################################################################
  # 1. ISME and MSE Table for Point Estimates
  #############################################################################
  
  # Define function to calculate the integrated square mean square error (ISME)
  # and mean square error (MSE)
  ISME_MSE <- function(estimate,  # Point estimate (vector, T x 1)
                       truth      # True value (vector, T x 1)
                       ){
    sum((estimate - truth)^2)/sum(truth^2)
  }
  
  # Define function to calculate the ISME for the eigenfunctions so that the direction
  # of the point estimate is aligned with the truth eigenfunction
  ISME_MSE_flip <- function(estimate,  # Point estimate (vector, T x 1)
                            truth      # True value (vector, T x 1)
                            ){
    if(sum((estimate - truth)^2)/sum(truth^2) > 
       sum(((-1 * estimate) - truth)^2)/sum(truth^2)){
      estimate <- -1 * estimate
    }
    sum((estimate - truth)^2)/sum(truth^2)
  }
  
  # Function to calculate the ISME/MSE values of the point estimates and generate table 
  ISME_MSE_Table <- function(){
    
    
    # a. ISME Mean Function
    
    ## i. Estimated mean function
    ISME_mu.hat <- mean(unlist(lapply(1:length(postSumms), function(run){
      ISME_MSE(postSumms[[run]]$mu$point_ests$mu.hat, truth_funcs$mu)
    })))
    
    ## i. MBD median mean 
    ISME_mu.MBD <- mean(unlist(lapply(1:length(postSumms), function(run){
      ISME_MSE(postSumms[[run]]$mu$point_ests$mu.MBD, truth_funcs$mu)
    })))
    
    
    # b. ISME Eigenfunctions
    for(k in 1:K){
      
      ## i. Estimated eigenfunction
      assign(paste("ISME_psi_", k, ".hat", sep = ""), mean(unlist(lapply(1:length(postSumms), function(run){
        ISME_MSE_flip(postSumms[[run]]$psi[[k]]$point_ests$psi.hat, truth_funcs$psi[[k]])
      }))))
      
      ## i. Estimated eigenfunction via covariance
      assign(paste("ISME_psi_", k, ".tilde", sep = ""), mean(unlist(lapply(1:length(postSumms), function(run){
        ISME_MSE_flip(postSumms[[run]]$psi[[k]]$point_ests$psi.tilde, truth_funcs$psi[[k]])
      }))))
      
      ## iii. MBD median eigenfunction
      assign(paste("ISME_psi_", k, ".MBD", sep = ""), mean(unlist(lapply(1:length(postSumms), function(run){
        ISME_MSE_flip(postSumms[[run]]$psi[[k]]$point_ests$psi.MBD, truth_funcs$psi[[k]])
      }))))
      
      ## iv. MVD median eigenfunction
      assign(paste("ISME_psi_", k, ".MVD", sep = ""), mean(unlist(lapply(1:length(postSumms), function(run){
        ISME_MSE_flip(postSumms[[run]]$psi[[k]]$point_ests$psi.MVD, truth_funcs$psi[[k]])
      }))))
    }
    
    
    # c. MSE Eigenvalues 
    for(k in 1:K){
      
      # i. Estimated eigenvalue
      assign(paste("ISME_rho_", k, ".hat", sep = ""), mean(unlist(lapply(1:length(postSumms), function(run){
        ISME_MSE(postSumms[[run]]$rho[[k]]$point_ests$rho.hat, truth_funcs$rho[[k]])
      }))))
      
      ## ii. Estimated eigenvalue via covariance 
      assign(paste("ISME_rho_", k, ".tilde", sep = ""), mean(unlist(lapply(1:length(postSumms), function(run){
        ISME_MSE(postSumms[[run]]$rho[[k]]$point_ests$rho.tilde, truth_funcs$rho[[k]])
      }))))
    }
    
    
    # d. Format data for table
    
    # Create vector of ISME and MSE values
    ISME_MSE_values <- c(ISME_mu.hat, ISME_mu.MBD)
    for(k in 1:K){
      ISME_MSE_values <- c(ISME_MSE_values, get(paste0("ISME_psi_", k, ".hat")),
                           get(paste0("ISME_psi_", k, ".tilde")),
                           get(paste0("ISME_psi_", k, ".MBD")),
                           get(paste0("ISME_psi_", k, ".MVD")))
    }
    for(k in 1:K){
      ISME_MSE_values <- c(ISME_MSE_values, get(paste0("ISME_psi_", k, ".hat")),
                           get(paste0("ISME_psi_", k, ".tilde")))
    }
    
    # Create vector of point estimate notation
    est_notation <- c("$\\widehat{\\mu}(t)$", "$\\widehat{m}\\{\\mu(t)\\}$")
    for(k in 1:K){
      est_notation <- c(est_notation, paste0("$\\widehat{\\psi}_", k, "(t)$"),
                        paste0("$\\widetilde{\\psi}_", k, "(t)$"),
                        paste0("$\\widehat{m}\\{\\psi_", k, "(t)\\}$"),
                        paste0("$\\widetilde{m}\\{\\psi_", k, "(t)\\}$"))
    }
    for(k in 1:K){
      est_notation <- c(est_notation, paste0("$\\widehat{\\rho}_", k, "$"),
                        paste0("$\\widetilde{\\rho}_", k, "$"))
    }
    
    # Create data.frame 
    ISME_MSE_values <- data.frame(est_notation = est_notation, value = ISME_MSE_values)
    colnames(ISME_MSE_values) <- c("Point Estimate", "")
    
    # Generate table 
    kbl(ISME_MSE_values, digits = 4, align = c('l', 'c')) %>%
      kable_classic(full_width = F) %>%
      group_rows("ISME", 1, 2 + K * 4, bold = T, underline = TRUE,
                 label_row_css = "text-align: right;", indent = F) %>%
      group_rows("MSE", (2 + K * 4) + 1, ((2 + K * 4) + 1) + K * 2 - 1, bold = T, underline = TRUE,
                 label_row_css = "text-align: right;", indent = F)
  }
  
  #############################################################################
  # 2. AR and NAR Table for Credible Intervals
  #############################################################################
  
  # Define function to calculate AR of credible interval with respect to the posterior 
  # estimates 
  AR <- function(data,         # Posterior estimates (matrix, M x T)
                 credible_int  # Credible interval estimates (list, length 2)
                                  # lower: lower bound (vector, T x 1)
                                  # upper: upper bound (vector, T x 1)
                 ){
    data.upper <- apply(data, 2, max)  # upper bounds of estimates
    data.lower <- apply(data, 2, min)  # lower bounds of estimates
    
    # Area of credible interval that lies inside the estimates 
    AR.lower <- ifelse(data.lower < credible_int$lower, 
                       credible_int$lower, data.lower)
    AR.upper <- ifelse(data.upper < credible_int$upper, 
                       data.upper, credible_int$upper)
    
    sum(AR.upper - AR.lower)/sum(data.upper - data.lower)
  }
  
  # Define function to calculate NAR of credible interval with respect to the posterior 
  # estimates 
  NAR <- function(data,         # Posterior estimates (matrix, M x T) 
                  credible_int  # Credible interval estimates (list, length 2)
                                  # lower: lower bound (vector, T x 1)
                                  # upper: upper bound (vector, T x 1)
                  ){
    data.upper <- apply(data, 2, max) # upper bounds of estimates
    data.lower <- apply(data, 2, min) # lower bounds of estimates
    
    negative_area <- sum(ifelse(  # Area of credible interval that lies outside the estimates
      data.upper < credible_int$upper, credible_int$upper - data.upper, 0)) +
      sum(ifelse(
        credible_int$lower < data.lower, data.lower - credible_int$lower, 0))
    
    negative_area/sum(data.upper - data.lower)
  }
  
  # Function to calculate the AR/NAR values of the credible intervals and generate Table 
  AR_NAR_Table <- function(){
    
    
    # a. AR Pointwise Parametric Credible Interval
    # Mean function
    AR_PP <- c(mean(
      unlist(lapply(1:length(postEsts), function(run){
        AR(postEsts[[run]]$model_components$mu.m, 
           postSumms[[run]]$mu$credible_ints$pointwise_parametric)
    }))))
    # Eigenfunctions
    for(k in 1:K){
      AR_PP <- c(AR_PP, mean(
        unlist(lapply(1:length(postEsts), function(run){
          AR(postEsts[[run]]$model_components$psi.m[[k]], 
             postSumms[[run]]$psi[[k]]$credible_ints$pointwise_parametric)
        }))))
    }
    AR_PP <- sprintf("%.3f", round(AR_PP, digits = 3))
    
    
    # b. AR Simultaneous Parametric Credible Interval
    # Mean function
    AR_SP <- c(mean(
      unlist(lapply(1:length(postEsts), function(run){
        AR(postEsts[[run]]$model_components$mu.m, 
           postSumms[[run]]$mu$credible_ints$simultaneous_parametric)
      }))))
    # Eigenfunctions
    for(k in 1:K){
      AR_SP <- c(AR_SP, mean(
        unlist(lapply(1:length(postEsts), function(run){
          AR(postEsts[[run]]$model_components$psi.m[[k]], 
             postSumms[[run]]$psi[[k]]$credible_ints$simultaneous_parametric)
        }))))
    }
    AR_SP <- sprintf("%.3f", round(AR_SP, digits = 3))
    
    
    # c. AR Pointwise Quantile Credible Interval
    # Mean function
    AR_PQ <- c(mean(
      unlist(lapply(1:length(postEsts), function(run){
        AR(postEsts[[run]]$model_components$mu.m, 
           postSumms[[run]]$mu$credible_ints$pointwise_quantile)
      }))))
    # Eigenfunctions
    for(k in 1:K){
      AR_PQ <- c(AR_PQ, mean(
        unlist(lapply(1:length(postEsts), function(run){
          AR(postEsts[[run]]$model_components$psi.m[[k]], 
             postSumms[[run]]$psi[[k]]$credible_ints$pointwise_quantile)
        }))))
    }
    AR_PQ <- sprintf("%.3f", round(AR_PQ, digits = 3))
    
    
    # d. AR Simultaneous Quantile Credible Interval'
    # Mean function
    AR_SQ <- c(mean(
      unlist(lapply(1:length(postEsts), function(run){
        AR(postEsts[[run]]$model_components$mu.m, 
           postSumms[[run]]$mu$credible_ints$simultaneous_quantile)
      }))))
    # Eigenfunctions
    for(k in 1:K){
      AR_SQ <- c(AR_SQ, mean(
        unlist(lapply(1:length(postEsts), function(run){
          AR(postEsts[[run]]$model_components$psi.m[[k]], 
             postSumms[[run]]$psi[[k]]$credible_ints$simultaneous_quantile)
        }))))
    }
    
    AR_SQ <- sprintf("%.3f", round(AR_SQ, digits = 3))
    
    
    # e. AR MBD Credible Interval
    # Mean function
    AR_MBD <- c(mean(
      unlist(lapply(1:length(postEsts), function(run){
        AR(postEsts[[run]]$model_components$mu.m, 
           postSumms[[run]]$mu$credible_ints$MBD)
      }))))
    # Eigenfunctions 
    for(k in 1:K){
      AR_MBD <- c(AR_MBD, mean(
        unlist(lapply(1:length(postEsts), function(run){
          AR(postEsts[[run]]$model_components$psi.m[[k]], 
             postSumms[[run]]$psi[[k]]$credible_ints$MBD)
        }))))
    }
    AR_MBD <- sprintf("%.3f", round(AR_MBD, digits = 3))
    
    
    # f. AR MVD Credible Interval
    AR_MVD <- c("")
    # Eigenfunctions
    for(k in 1:K){
      AR_MVD <- c(AR_MVD, sprintf("%.3f", round(mean(
        unlist(lapply(1:length(postEsts), function(run){
          AR(postEsts[[run]]$model_components$psi.m[[k]], 
             postSumms[[run]]$psi[[k]]$credible_ints$MVD)
        }))), 3)))
    }
    
    
    # g. NAR Pointwise Parametric Credible Interval
    # Mean function
    NAR_PP <- c(mean(
      unlist(lapply(1:length(postEsts), function(run){
        NAR(postEsts[[run]]$model_components$mu.m, 
           postSumms[[run]]$mu$credible_ints$pointwise_parametric)
      }))))
    # Eigenfunctions
    for(k in 1:K){
      NAR_PP <- c(NAR_PP, mean(
        unlist(lapply(1:length(postEsts), function(run){
          NAR(postEsts[[run]]$model_components$psi.m[[k]], 
             postSumms[[run]]$psi[[k]]$credible_ints$pointwise_parametric)
        }))))
    }
    NAR_PP <- sprintf("%.3f", round(NAR_PP, digits = 3))
    
    
    # h. NAR Simultaneous Parametric Credible Interval
    # Mean function
    NAR_SP <- c(mean(
      unlist(lapply(1:length(postEsts), function(run){
        NAR(postEsts[[run]]$model_components$mu.m, 
           postSumms[[run]]$mu$credible_ints$simultaneous_parametric)
      }))))
    # Eigenfunctions
    for(k in 1:K){
      NAR_SP <- c(NAR_SP, mean(
        unlist(lapply(1:length(postEsts), function(run){
          NAR(postEsts[[run]]$model_components$psi.m[[k]], 
             postSumms[[run]]$psi[[k]]$credible_ints$simultaneous_parametric)
        }))))
    }
    NAR_SP <- sprintf("%.3f", round(NAR_SP, digits = 3))
    
    
    # i. NAR Simultaneous Quantile Credible Interval
    # Mean function
    NAR_SQ <- c(mean(
      unlist(lapply(1:length(postEsts), function(run){
        NAR(postEsts[[run]]$model_components$mu.m, 
           postSumms[[run]]$mu$credible_ints$simultaneous_quantile)
      }))))
    # Eigenfunctions
    for(k in 1:K){
      NAR_SQ <- c(NAR_SQ, mean(
        unlist(lapply(1:length(postEsts), function(run){
          NAR(postEsts[[run]]$model_components$psi.m[[k]], 
             postSumms[[run]]$psi[[k]]$credible_ints$simultaneous_quantile)
        }))))
    }
    NAR_SQ <- sprintf("%.3f", round(NAR_SQ, digits = 3))
    
    
    # j. Format data.frame for table generation
    func_names <- c("$\\mu(t)$")
    for(k in 1:K){
      func_names <- c(func_names, paste0("$\\psi_", k, "(t)$"))
    }
    
    AR_NAR_value <- data.frame(cbind(func_names, AR_PP, AR_SP, AR_PQ, AR_SQ, 
                                     AR_MBD, AR_MVD, NAR_PP, NAR_SP, NAR_SQ))
    AR_NAR_value[AR_NAR_value == 0] <- "0.000"
    
    colnames(AR_NAR_value) <- c("$g(t)$", "$P^p_{1 - \\alpha}\\{g(t)\\}$",
                                "$P^s_{1 - \\alpha}\\{g(t)\\}$",
                                "$Q^p_{1 - \\alpha}\\{g(t)\\}$", 
                                "$Q^s_{1 - \\alpha}\\{g(t)\\}$", 
                                "$D_{1 - \\alpha}\\{g(t)\\}$",
                                "$D^{\\star}_{1 - \\alpha}\\{g(t)\\}$", 
                                "$P^p_{1 - \\alpha}\\{g(t)\\}$", 
                                "$P^s_{1 - \\alpha}\\{g(t)\\}$",
                                "$Q^s_{1 - \\alpha}\\{g(t)\\}$")
    
    # Generate table
    kbl(AR_NAR_value, align = c('l', rep('c', 9))) %>%
      kable_classic(full_width = F) %>%
      add_header_above(c(" " = 1, "AR" = 6, "NAR" = 3), bold = T)
  }
  
  ##### Plot table #####
  
  if(tableIndex == 1){
    ISME_MSE_Table()
  } else {
    AR_NAR_Table()
  }
}