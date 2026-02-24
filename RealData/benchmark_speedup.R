# New File: benchmark_speedup.R (Reproduces Table 1-style speedup)

library(tidyverse)
library(microbenchmark)

source("bootstrap_functions.R")

# Benchmark function
benchmark_bootstrap <- function(n_sizes = c(1e3, 1e4, 1e5, 1e6), B = 10000) {
  results <- data.frame()
  
  for (n in n_sizes) {
    x <- rnorm(n)
    
    # Serial
    time_serial <- microbenchmark(
      studentized_bootstrap(x, B = B, parallel = FALSE),
      times = 3
    )$time / 1e9  # seconds
    
    # Parallel
    time_parallel <- microbenchmark(
      studentized_bootstrap(x, B = B, parallel = TRUE),
      times = 3
    )$time / 1e9
    
    speedup <- mean(time_serial) / mean(time_parallel)
    
    results <- rbind(results, data.frame(
      n = n,
      serial_sec = round(mean(time_serial), 2),
      parallel_sec = round(mean(time_parallel), 2),
      speedup = round(speedup, 1)
    ))
  }
  
  print(results)
  results
}

cat("Running benchmarks (may take several minutes for n=10^6)...\n")
benchmark_bootstrap()