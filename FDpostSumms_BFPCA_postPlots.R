postPlots <- function(data,          # Sample of posterior estimates, g^{(m)}(t) (matrix, M x T)
                      k,             # Index for the eigenfunctions taking values 1,...K, but when 
                                     # plotting the mean function takes value 0 (integer)
                      t,             # Functional domain grid (vector, T x 1)
                      postSumms,     # list returned from postSumms() function containing the estimated 
                                     # posterior summaries for the model components
                      point_est,     # Measure of central tendency used for plot. Takes values 
                                     # c("hat", "tilde", "MBD", "MVD")
                      credible_int,  # Measure of spread used for plot. Takes values
                                     # c("parametric", "quantile", "depth")
                      title,         # Title given to the plot (character)
                      ylab,          # Y-axis label (character)
                      xlab           # X-axis label (character)
                      ){
  #############################################################################
  ## Description: Function for generating the individual plots of the functional 
  ##              posterior estimates calculated for a given simulation run displayed in 
  ##              Figures 3 and 4 of "Functional Depth Posterior Summaries for Bayesian 
  ##              Functional Principal Component Analysis" by Boland et al. (2022). 
  ##              The user can formally select which point estimate (solid black line) 
  ##              and credible intervals (shaded area and dashed black lines) they 
  ##              wish to display using the arguments point_est and credible_int, 
  ##              respectively. To note, the MVD credible interval cannot be calculated
  ##              for the mean function, and as a result, the solid dashed lines 
  ##              are not displayed when k = 0 and credible_int = "depth". 
  ## Args:        (see above)
  ## Returns:     list()
  ## postPlots Outline:
  ##              1. Clean and format data
  ##              2. Plot the visualization
  #############################################################################
  
  # Install missing packages
  list.of.packages <- c("tidyverse", "reshape2")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.packages) 
  
  # Load packages
  library(tidyverse)
  library(reshape2)
  
  #############################################################################
  # 1. Clean and format data
  #############################################################################
  
  # Format the sample posterior estimates from wide to long format
  data <- data.frame(t(data))
  data$Time <- t
  data <- data %>%
    melt(id.vars = c("Time"))
  
  # Obtain set of point estimates 
  if(k == 0){
    postSumms <- postSumms$mu
  } else {
    postSumms <- postSumms$psi[[k]]
  }
  
  # Obtain the point estimate for plotting 
  if(point_est == "hat" && k != 0){
    point_estimate <- postSumms$point_ests$psi.hat 
  } else if (point_est == "tilde" && k != 0) {
    point_estimate <- postSumms$point_ests$psi.tilde 
  } else if (point_est == "MBD" && k != 0) {
    point_estimate <- postSumms$point_ests$psi.MBD 
  } else if (point_est == "MVD" && k != 0) {
    point_estimate <- postSumms$point_ests$psi.MVD 
  } else if (point_est == "hat" && k == 0) {
    point_estimate <- postSumms$point_ests$mu.hat 
  } else if (point_est == "MBD" && k == 0) {
    point_estimate <- postSumms$point_ests$mu.MBD 
  }
  
  # Format data for point estimates
  point_estimate <- data.frame(Time = t, point_estimate = point_estimate)
  
  # Obtain upper and lower bounds of credible intervals for plotting
  if(credible_int == "parametric"){
    pointwise.lower <- postSumms$credible_ints$pointwise_parametric$lower
    pointwise.upper <- postSumms$credible_ints$pointwise_parametric$upper
    
    simultaneous.lower <- postSumms$credible_ints$simultaneous_parametric$lower
    simultaneous.upper <- postSumms$credible_ints$simultaneous_parametric$upper
    
    fill.color <- "#00AFBB"
  } else if (credible_int == "quantile") {
    pointwise.lower <- postSumms$credible_ints$pointwise_quantile$lower
    pointwise.upper <- postSumms$credible_ints$pointwise_quantile$upper
    
    simultaneous.lower <- postSumms$credible_ints$simultaneous_quantile$lower
    simultaneous.upper <- postSumms$credible_ints$simultaneous_quantile$upper
    
    fill.color <- "#FC4E07"
  } else if (credible_int == "depth" && k != 0) {
    pointwise.lower <- postSumms$credible_ints$MBD$lower
    pointwise.upper <- postSumms$credible_ints$MBD$upper
    
    simultaneous.lower <- postSumms$credible_ints$MVD$lower
    simultaneous.upper <- postSumms$credible_ints$MVD$upper
    
    fill.color <- "#E7B800"
  } else if (credible_int == "depth" && k == 0) {
    pointwise.lower <- postSumms$credible_ints$MBD$lower
    pointwise.upper <- postSumms$credible_ints$MBD$upper
    
    fill.color <- "#E7B800"
  } 
  
  # Format data for credible intervals 
  if(credible_int == "depth" && k == 0){
    point.ci <- data.frame(cbind(Time = t, pointwise.lower = pointwise.lower,
                                 pointwise.upper = pointwise.upper))
  } else {
    point.ci <- data.frame(cbind(Time = t, pointwise.lower = pointwise.lower,
                                 pointwise.upper = pointwise.upper))
    
    simlt.ci <- data.frame(cbind(Time = t, simultaneous.lower = simultaneous.lower, 
                                 simultaneous.upper = simultaneous.upper))
    simlt.ci <- simlt.ci %>%
      melt(id.vars = c("Time"))
  }
  
  #############################################################################
  # 2. Plot the visualization
  #############################################################################
  
  if(credible_int == "depth" && k == 0){  # Plot functional depth intervals for mean function (omits MVD interval)
    ggplot() +
      geom_line(data = data,  # plot posterior estimates 
                mapping = aes(x = Time, y = value, color = variable)) +
      geom_ribbon(data = point.ci, mapping = aes(x = Time, ymin = pointwise.lower,  # plot pointwise or MBD intervals
                                                 ymax = pointwise.upper),
                  fill = fill.color, alpha = 0.5) + 
      theme_bw() + 
      labs(title = title, y = ylab, x = "Time (t)") +
      theme(legend.position = "none", 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 12)) +
      scale_color_manual(values = c(rep("gray88", 4000))) +  # Color of the posterior estimates (grey lines)
      scale_x_continuous(expand = c(0, 0)) +
      geom_line(data = point_estimate, mapping = aes(x = Time, y = point_estimate), color = "black", # plot point estimate
                size = 0.9)
  } else { # All other plots 
    ggplot() +
      geom_line(data = data,  # plot posterior estimates 
                mapping = aes(x = Time, y = value, color = variable)) +
      geom_ribbon(data = point.ci, mapping = aes(x = Time, ymin = pointwise.lower,  # plot pointwise or MBD intervals
                                                 ymax = pointwise.upper),
                  fill = fill.color, alpha = 0.5) + 
      theme_bw() + 
      labs(title = title, y = ylab, x = "Time (t)") +
      theme(legend.position = "none", 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 12)) +
      scale_color_manual(values = c(rep("gray88", 4000))) +  # Color of the posterior estimates (grey lines)
      scale_x_continuous(expand = c(0, 0)) +
      geom_line(data = point_estimate, mapping = aes(x = Time, y = point_estimate), color = "black", # plot point estimate
                size = 0.9) +
      geom_line(data = simlt.ci, # plot simultaneous or MVD intervals 
                mapping = aes(x = Time, y = value, linetype = variable), 
                size = 0.9) +
      scale_linetype_manual(values = c(2, 2))
  }
}

