
.onLoad <- function(libname, pkgname) {
}


#' Install BRACoD in python
#'
#' Uses pip to install the latest BRACoD release in python. You might need
#' to specify a python environment with either reticulate::use_virtualenv or
#' reticulate::use_condaenv. 
#' @param method passed to reticulate::py_install
#' @param conda passed to reticulate::py_install
#' @return no return value
#' @export
install_bracod <- function(method = "auto", conda = "auto") {
  reticulate::py_install("BRACoD", method = method, conda = conda, pip=TRUE)
}


#' Simulate microbiome counts
#' 
#' Each bacteria's absolute abundance is simulated from a lognormal distribution.
#' Then, convert each sample to relative abundance, and simulate sequencing counts
#' using a multinomial distribution, based on the desired number of reads and the 
#' simulated relative abundances. This also simulates an environmental variable that
#' is produced by some of the bacteria.
#' @param df A dataframe of OTU counts that is a model for data simulation. Samples are rows and bacteria are columns.
#' @param n_contributors the number of bacteria that are to contribute to your environmental variable.
#' @param coeff_contributor the average of the distribution used to simulate the contribution coefficient.
#' @param min_ab_contributor The minimum log relative abundance, averaged across samples, to include a bacteria
#' @param sd_Y the standard deviation of the simulated environmental variable
#' @param n_reads the number of reads to be simulated per sample
#' @param var_contributor If you use a uniform distribution, this is the range of the distribution, with a normal distribution it is the variance used to simulate the contribution coefficient.
#' @param use_uniform use a uniform distribution to simulate the contribution coefficient. Alternative is the normal distribution.
#' @param n_samples_use number of microbiome samples to simulate. If NULL, uses the same number of samples as in your dataframe
#' @param corr_value the bacteria-bacteria correlation value you want to include in the simulation
#' @param return_absolute returns the abosulte abundance values instead of the simulated microbiome counts
#' @param seed random seed for reproducibility
#' @return a list containing 1) the simulated count data 2) the simulated environmental variable and 3) the simulated contribution coefficients
#' @export
simulate_microbiome_counts <- function(df, n_contributors = 20, coeff_contributor = 0.0, min_ab_contributor = -9, sd_Y = 1.0, n_reads = 100000, var_contributor = 5.0, use_uniform = TRUE, n_samples_use = NULL, corr_value = NULL,  return_absolute = FALSE, seed = NULL) {
  BRACoD <- reticulate::import("BRACoD")
  results <- BRACoD$simulate_microbiome_counts(df, n_contributors = n_contributors, coeff_contributor = coeff_contributor, min_ab_contributor = min_ab_contributor, sd_Y = sd_Y, n_reads = n_reads, var_contributor = var_contributor, use_uniform = use_uniform, n_samples_use = n_samples_use, corr_value = corr_value, return_absolute = return_absolute, seed = seed)
  return(list(sim_counts=results[[1]],sim_y=results[[2]],contributions=results[[3]]))
}


#' Normalize OTU counts and add a pseudo count
#'
#' BRACoD requires relative abundance and cannot handle zeros, so this function
#' adds a small pseudo count (1/10th the smallest non-zero value).
#' @param df_counts A dataframe of OTU counts. Samples are rows and bacteria are columns.
#' @return a dataframe of relative abundance data
#' @export
scale_counts <- function(df_counts) {
  BRACoD <- reticulate::import("BRACoD")
  # Check if this is legitimately counts data
  stopifnot(all(apply(df_counts, 1, function(x) all(x == as.integer(x)))))
  
  # Frequently, R and python conversion results in a "double" dataframe that has counts data. We need to fix that
  df_counts <- t(apply(df_counts, 1, function(x) as.integer(x)))

  return(BRACoD$scale_counts(df_counts))
}

#' Threshold your microbiome counts data
#'
#' This function removes samples below a minimum counts and bacteria below a minimum log abundance.
#' Run this before running BRACoD because the algorithm does not perform well when there are many
#' low abundance bacteria that are only present in a few samples.
#' @param df_counts A dataframe of OTU counts. Samples are rows and bacteria are columns.
#' @param min_counts threshold samples with fewer than this many counts
#' @param min_ab threshold bacteria whose average log abundance is below this
#' @return a dataframe of microbiome counts
#' @export
threshold_count_data <- function(df_counts, min_counts = 1000, min_ab=1e-4) {
  BRACoD <- reticulate::import("BRACoD")

  # Check if this is legitimately counts data
  stopifnot(all(apply(df_counts, 1, function(x) all(x == as.integer(x)))))

  # Frequently, R and python conversion results in a "double" dataframe that has counts data. We need to fix that
  df_counts_c <- t(apply(df_counts, 1, function(x) as.integer(x)))
  df_counts_c <- data.frame(df_counts_c)
  rownames(df_counts_c) <- rownames(df_counts)
  colnames(df_counts_c) <- colnames(df_counts)

  return(BRACoD$threshold_count_data(df_counts, min_counts=min_counts, min_ab=min_ab))
}


#' Run the main BRACoD algorithm
#'
#' Uses pymc3 to sample the posterior of the model to determine bacteria that are
#' associated with your environmental variable.
#' @param df_relab A dataframe of relative microbiome abundances. Samples are rows and bacteria are columns.
#' @param env_var the environmental variable you are evaluating. You need 1 measurement associated with each sample.
#' @param n_sample number of posterior samples.
#' @param n_burn number of burn-in steps before actual sampling stops.
#' @param njobs number of parallel MCMC chains to run.
#' @return the pymc trace object which holds the samples of the posterior distribution
#' @export
#' @examples
#' \dontrun{
#' data(obesity)
#' r <- simulate_microbiome_counts(obesity)
#' sim_counts <- r[[1]]
#' sim_y <- r[[2]]
#' contributions <- r[[3]]
#' sim_relab <- scale_counts(sim_counts)
#' trace <- run_bracod(sim_relab, sim_y, n_sample = 1000, n_burn=1000, njobs=4)
#' }
run_bracod <- function(df_relab, env_var, n_sample=1000, n_burn=1000, njobs=4) {
  BRACoD <- reticulate::import("BRACoD")
  return(BRACoD$run_bracod(df_relab, env_var, n_sample = n_sample, n_burn=n_burn, njobs=njobs))
}


#' Summarize the results of BRACoD
#'
#' This summarizes the trace object that run_bracod() returns. It returns a dataframe
#' that contains two parameters of interest, the average inclusion (p) and the average
#' coefficient (beta), telling you the association between that bacteria and the environmental
#' variable
#' @param trace the pymc3 object that is the output of run_bracod()
#' @param taxon_names optional, a list of names of the bacteria to include in the results
#' @param cutoff this is the cutoff on the average inclusion for inclusion. We reccomend a value of 0.3, but you can lower the value to include less confident taxon or raise the cutoff to exclude them.
#' @return a dataframe with information about the bacteria that BRACoD identified
#' @export
#' @examples
#' \dontrun{
#' trace <- run_bracod(sim_relab, sim_y, n_sample = 1000, n_burn=1000, njobs=4)
#' df_summary <- summarize_trace(trace, colnames(sim_relab))
#' }
summarize_trace <- function(trace, taxon_names=NULL, cutoff=0.3) {
  BRACoD <- reticulate::import("BRACoD")
  df <- BRACoD$summarize_trace(trace, taxon_names, cutoff)
  # Python to R numbering
  df$taxon_num <- df$taxon_num + 1
  return(df)
}

#' Score the results of BRACoD
#'
#' This calculate the precision, recall and F1 of your BRACoD results if you know
#' the ground truth, ie. if this is simulated data.
#' @param taxon_identified a list of integers corresponding to the indicies of the taxon you identified with BRACoD
#' @param taxon_actual a list of integers corresponding to the indicies of the taxon that truely contribute to butyrate levels
#' @return a list containing 1) the precision 2) the recall 3) the f1 metric
#' @export
#' @examples
#' \dontrun{
#' df_summary <- summarize_trace(trace, colnames(sim_relab))
#' taxon_identified <- df_summary$taxon
#' taxon_actual <- which(contributions != 0)
#' 
#' r <- score(taxon_identified, taxon_actual)
#' 
#' precision <- r[[1]]
#' recall <- r[[2]]
#' f1 <- r[[3]]
#' 
#' print(sprintf("Precision: %.2f, Recall: %.2f, F1: %.2f",precision, recall, f1))
#' }
score <- function(taxon_identified, taxon_actual) {
  BRACoD <- reticulate::import("BRACoD")
  results <- BRACoD$score(taxon_identified, taxon_actual)
  return(list(precision=results[[1]],recall=results[[2]],f1=results[[3]]))
}


#' Remove NULL values in your OTU and environmental variable
#'
#' This will remove samples that are NULL in the environmental variable, as well as
#' the corresponding samples in your relative abundance data.
#' @param df_relab microbiome relative abundance data in a dataframe
#' @param Y values of the environmental variable
#' @return a list containing 1) the relative abundance data and 2) the Y values
#' @export
remove_null <- function(df_relab, Y) {
  BRACoD <- reticulate::import("BRACoD")
  results <- BRACoD$remove_null(df_relab, Y)
  return(list(df_rel=results[[1]],Y=results[[2]]))
}


#' Perform convergence tests on the p and beta variables
#'
#' You may get errors are divergence of some variables after pymc3 samples the posterior.
#' We are not overly concerned about some of the variables, such as the variance, rather
#' we are really interested in the inclusion probabilities (p) and contribution coefficients
#' (beta). The convergence tests that are included here focus on evaluating those two variables.
#' @param trace the output of run_bracod()
#' @param df_relab the microbiome relative abundance
#' @return no return value
#' @export
convergence_tests <- function(trace, df_relab) {
  BRACoD <- reticulate::import("BRACoD")
  BRACoD$convergence_tests(trace, df_relab)
}

