# Compile and load C++ likelihood functions
# Run this script once to compile the C++ code

# Check for required packages
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  install.packages("Rcpp")
}
if (!requireNamespace("RcppArmadillo", quietly = TRUE)) {
  install.packages("RcppArmadillo")
}

library(Rcpp)
library(RcppArmadillo)

# Compile the C++ source file
# This creates the functions: NIMR_cpp, lglik_cpp, gradLLS_cpp, hessLLS_cpp
sourceCpp("src/likelihood_functions.cpp")

cat("C++ functions compiled successfully!\n")
cat("Available functions: NIMR_cpp, lglik_cpp, gradLLS_cpp, hessLLS_cpp\n")
