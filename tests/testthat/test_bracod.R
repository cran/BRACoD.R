skip_if_no_bracod <- function() {
  have_bracod <- reticulate::py_module_available("BRACoD")
  if (!have_bracod)
    skip("bracod not available for testing")
}

test_that("remove null works", {
   data(obesity)
   skip_if_no_bracod()
   r <- simulate_microbiome_counts(df_counts_obesity)

   sim_counts <- r[[1]]
   sim_y <- r[[2]]
   sim_relab <- scale_counts(sim_counts)
   r <- remove_null(sim_relab, sim_y)
   
   expect_equal(nrow(r[[1]]), length(r[[2]]))
})

test_that("bracod runs", {
    data(obesity)
    skip_if_no_bracod()
    r <- simulate_microbiome_counts(df_counts_obesity)

    sim_counts <- r[[1]]
    sim_y <- r[[2]]
    contributions <- r[[3]]
    sim_relab <- scale_counts(sim_counts)
    trace <- run_bracod(sim_relab, sim_y, n_sample = 1000, n_burn=1000, njobs=4)
})
