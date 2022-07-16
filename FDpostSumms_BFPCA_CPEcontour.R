CPEcontour <- function(data,          # Sample of posterior estimates, g^{(m)}(t) (matrix, M x T)
                        Cov,          # Sample of posterior covariance estimates (list, length M)
                        CPE,          # Type of CPE and point estimate to be used in plot. 
                                      # Takes values c("MBD", "MVD"),
                        alpha_cont,   # Vector of alpha-levels to be included in the contour between 
                                      # values c(0, 1)
                        t,            # Functional domain grid (vector, T x 1)
                        title,        # Title given to the plot (character)
                        ylab,         # Y-axis label (character)
                        xlab          # X-axis label (character)
){
  #############################################################################
  ## Description: Function for generating the CPE contour plots of the functional 
  ##              posterior estimates calculated for a given simulation run displayed in 
  ##              Figures 3, 4 and S3 of "Central Posterior Envelopes for Bayesian 
  ##              Functional Principal Component Analysis" by Boland et al. (2022). 
  ##              The user can formally select which CPE type (MVD or MBD) and 
  ##              respective  MBD/MVD median (solid black line) and CPEs (shaded areas) 
  ##              they wish to display at different alpha-levels using the arguments 
  ##              CPE and alpha_cont, respectively. To note, the MVD-CPE cannot be calculated
  ##              for the mean function.
  ## Args:        (see above)
  ## Returns:     list()
  ## CPEcontour Outline:
  ##              1. Clean and format data
  ##              3. Obtain legend
  ##              3. Plot the visualization
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("tidyverse", "reshape2", "viridis", "cowplot",
                        "latex2exp")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.packages) 
  
  # Load packages
  library(tidyverse)
  library(reshape2)
  library(viridis)
  library(cowplot)
  library(latex2exp)
  
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
  
  #############################################################################
  # 1. Clean and format data
  #############################################################################
  
  # Order the posterior sample based on MBD or MVD
  if(CPE == "MBD"){
    ranking <- MBD(t(data))
  } else if (CPE == "MVD"){
    ranking <- MVD(Cov)
  }
  ranking <- ranking[order(-ranking$rank),]
  ordering_index <- ranking$id
  data_ordered <- as.matrix(data[ordering_index, ])
  
  # Obtain point estimate data
  point_est <- data.frame(Time = t, est = data_ordered[1, ]) 
  
  # Calculate the CPE regions 
  total_alpha <- length(alpha_cont)
  CPE_regions <- lapply(X = 1:total_alpha, function(X){
     
    # Calculate upper bounds
    upper <- apply(data_ordered[1:(floor(alpha_cont[X] * nrow(data_ordered))),], 2, max)
    
    # Calculate lower bounds
    lower <- apply(data_ordered[1:(floor(alpha_cont[X] * nrow(data_ordered))),], 2, min)
    
    return(list(upper = upper, lower = lower))
  })
  
  # Turn the data from wide to long format
  data <- data.frame(t(data))
  data$Time <- t
  data <- data %>%
    melt(id.vars = c("Time"))
  
  # Generate datasets to plot CPE contours 
  contour_min <- data.frame(cbind(Time = t, upper = CPE_regions[[1]]$upper,
                                  lower = CPE_regions[[1]]$lower))
  
  contour_bands <- lapply(2:total_alpha, function(X){
    contour_upper <- data.frame(cbind(Time = t, upper = CPE_regions[[X]]$upper,
                                      lower = CPE_regions[[X - 1]]$upper))
    
    contour_lower <- data.frame(cbind(Time = t, upper = CPE_regions[[X - 1]]$lower,
                                      lower = CPE_regions[[X]]$lower))
    
    return(list(contour_upper = contour_upper, contour_lower = contour_lower))
  })
  
  #############################################################################
  # 2. Obtain legend
  #############################################################################
  
  # Colors used in contour plot
  cols <- viridis::viridis(n = total_alpha)
  
  # Generate a dataset and plot to obtain legend
  df <- data.frame(x = rnorm(1000, sd = 100), y = rnorm(1000, sd = 100))
  grid <- ggplot(df, aes(x = x, y = y)) +
    geom_density_2d_filled(bins = 10) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title=element_text(size=12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 13),
          legend.background = element_rect(colour = 'black', fill = 'white', 
                                           linetype='solid'),
          axis.text=element_text(size=12)) +
    scale_fill_manual(values =  cols, 
                      labels = c(as.character(alpha_cont))) +
    guides(fill=guide_legend(title=TeX(r'($\alpha$ level)')))
  legend <- get_legend(grid)

  #############################################################################
  # 3. Plot the visualization
  #############################################################################
  
  contour_plot <- ggplot() +
    geom_line(data = data, 
              mapping = aes(x = Time, y = value, color = variable), alpha = 0.3) +
    scale_color_manual(values = c(rep("gray70", 4000))) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_bw() + 
    labs(title = title, y = ylab, x = xlab) +
    theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 12)) +
    geom_ribbon(data = contour_min, mapping = aes(x = Time, ymin = lower,
                                                 ymax = upper),
                fill = cols[1], alpha = 0.8)
  
  for(i in 2:total_alpha){
    contour_plot <- contour_plot + 
      geom_ribbon(data = contour_bands[[i - 1]]$contour_upper, 
                  mapping = aes(x = Time, ymin = lower, ymax = upper),
                  fill = cols[i], alpha = 0.8) +
      geom_ribbon(data = contour_bands[[i - 1]]$contour_lower, 
                  mapping = aes(x = Time, ymin = lower, ymax = upper),
                  fill = cols[i], alpha = 0.8) 
  }
  
  contour_plot <- contour_plot + 
    geom_line(data = point_est, mapping = aes(x = Time, y = est), color = "black",
              size = 0.9) 
  
  ggdraw() +
    draw_plot(contour_plot, x = 0, y = 0, width = 0.9 - 0.01, height = 1) +
    draw_grob(legend, x = 0.9, y = 0, width = 0.1 - 0.01, height = 1) 
}

