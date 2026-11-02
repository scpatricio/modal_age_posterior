`.sourceCpp_1_DLLInfo` <- dyn.load('/Users/silca/Desktop/fertility/PDF_of_the_mode/Code/R/random_walk/Results/.rcpp_cache/sourceCpp-aarch64-apple-darwin23-1.1.2/sourcecpp_113b47128f35c/sourceCpp_3.so')

rw_grid_cpp <- Rcpp:::sourceCppFunction(function(D, E, order, log10_lambdas, drop, max_iter = 50L, tol = 1e-7) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_rw_grid_cpp')
rw_draw_d_cpp <- Rcpp:::sourceCppFunction(function(eta, L, cell, order) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_rw_draw_d_cpp')
mode_stats_cpp <- Rcpp:::sourceCppFunction(function(d, ages, c0, c1, major) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_mode_stats_cpp')

rm(`.sourceCpp_1_DLLInfo`)
