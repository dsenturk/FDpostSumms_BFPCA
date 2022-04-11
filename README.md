CONTENTS OF THIS FOLDER ——————————————
* FDpostSumms_BFPCA_tutorial.Rmd : A step-by-step of implementation of the posterior estimation using the BFPCA model and subsequent calculation of traditional and proposed functional depth based posterior summaries described in "Functional Depth Posterior Summaries for Bayesian Functional Principal Component Analysis" by Boland et al. (2022).

* FDpostSumms_BFPCA_simulateData.R : Function for simulating functional data from a standard FPCA model for posterior estimation of BFPCA model components. 

* FDpostSumms_BFPCA_MCMC.R : Function for estimation of the BFPCA model posterior parameters using a Gibbs sampler and posterior model components including the mean function, covariance function, eigenvalues, and eigenfunctions. 

* FDpostSumms_BFPCA_postSumms.R : Function for calculating the traditional and proposed functional depth posterior summaries including point estimates and credible intervals from the estimates of the posterior model components returned from the MCMC() function. 

* FDpostSumms_BFPCA_postPlots.R : Function for visualization of the posterior estimates and summaries of the mean function and eigenfunctions. 

* FDpostSumms_BFPCA_simTable.R : Function for generating tables containing the performance measures for evaluation of the posterior summaries including point estimate evaluation using ISME/MSE and credible interval evaluation using AR/NAR.

INTRODUCTION ——————————————

The contents of this folder allow for implementation of the described in "Functional Depth Posterior Summaries for Bayesian Functional Principal Component Analysis" by et al. (2022). Users can simulate a sample data frame (FDpostSumms_BFPCA_simulateData.R) and apply the proposed BFPCA model posterior estimation (FDpostSumms_BFPCA_MCMC.R). Further, we include tools to calculate posterior summaries of the model components (FDpostSumms_BFPCA_postSumms.R), allowing users to describe the central tendency and spread of the functional posterior samples. Detailed instructions on how to perform the aforementioned procedures, visualize results (FDpostSumms_BFPCA_postPlots.R), and generate tables containing performance measures of the posterior summaries (FDpostSumms_BFPCA_simTable.R) are included in FDpostSumms_BFPCA_tutorial.Rmd.

REQUIREMENTS ——————————————

The included R programs require R  (R Core Team, 4.1.3) and the packages listed in the aformentioned scripts mentioned above.

INSTALLATION ——————————————

Load the R program files into the global environment and install required packages using steps detailed in FDpostSumms_BFPCA_tutorial.Rmd.
