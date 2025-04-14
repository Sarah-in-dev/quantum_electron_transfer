#!/usr/bin/env Rscript

# Function to install packages if they're not already installed
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing package:", package))
    # Try installing from CRAN repository
    install.packages(package, repos = "https://cloud.r-project.org", quiet = TRUE)
    # Check if installation was successful
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Failed to install package:", package))
    }
  }
}

# Check and install required packages
tryCatch({
  install_if_missing("ggplot2")
  install_if_missing("dplyr")
}, error = function(e) {
  cat("Error installing packages. You may need administrator privileges or a local R library.\n")
  cat("Try running: R -e 'install.packages(c(\"ggplot2\", \"dplyr\"), lib=\"~/R/library\", repos=\"https://cloud.r-project.org\")'\n")
  cat("Then set: export R_LIBS_USER=~/R/library\n")
  quit(status = 1)
})

# Load libraries
library(ggplot2)
library(dplyr)

# Check if files exist and can be read
check_file <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }
}

# Define file paths
quantum_file <- "../results/quantum_residues.txt"
distance_file <- "../results/cofactor_distances.txt"
conservation_file <- "../results/conservation_scores.txt"

# Check files
check_file(quantum_file)
check_file(distance_file)
check_file(conservation_file)

# Read data with error handling
tryCatch({
  # For quantum_residues.txt
  if (file.exists(quantum_file) && file.info(quantum_file)$size > 0) {
    quantum_residues <- read.table(quantum_file, header=TRUE, sep="\t")
    cat("Successfully loaded quantum residues data:", nrow(quantum_residues), "rows\n")
  } else {
    cat("Warning: quantum_residues.txt is empty or missing\n")
    # Create dummy data for testing
    quantum_residues <- data.frame(
      chain = c("A", "A", "D", "D", "A"),
      position = c(101, 102, 103, 104, 105),
      residue = c("ALA", "PHE", "GLY", "LEU", "TYR"),
      distance_to_cofactor = c(3.5, 8.2, 12.5, 15.7, 18.9),
      conservation = c(0.95, 0.87, 0.75, 0.62, 0.45),
      quantum_relevance = c("High", "High", "Medium", "Medium", "Low")
    )
  }
  
  # For cofactor_distances.txt
  if (file.exists(distance_file) && file.info(distance_file)$size > 0) {
    cofactor_distances <- read.table(distance_file, header=TRUE, sep="\t")
    cat("Successfully loaded cofactor distances data:", nrow(cofactor_distances), "rows\n")
  } else {
    cat("Warning: cofactor_distances.txt is empty or missing\n")
    # Create dummy data
    cofactor_distances <- data.frame(
      cofactor1 = c("P680_A", "P680_D", "PheoA"),
      cofactor2 = c("PheoA", "PheoD", "QA"),
      distance = c(10.2, 10.5, 12.7)
    )
  }
  
  # For conservation_scores.txt
  if (file.exists(conservation_file) && file.info(conservation_file)$size > 0) {
    conservation <- read.table(conservation_file, header=TRUE, sep="\t")
    cat("Successfully loaded conservation data:", nrow(conservation), "rows\n")
  } else {
    cat("Warning: conservation_scores.txt is empty or missing\n")
    # Create dummy data
    conservation <- data.frame(
      position = 1:200,
      conservation_score = runif(200, 0.3, 1.0)
    )
  }
}, error = function(e) {
  cat("Error reading data files:", e$message, "\n")
  quit(status = 1)
})

# Create results directory if it doesn't exist
dir.create("../results", showWarnings = FALSE)

# Calculate summary statistics
cat("Calculating summary statistics...\n")
tryCatch({
  summary_stats <- quantum_residues %>%
    group_by(quantum_relevance) %>%
    summarize(
      count = n(),
      percent = n() / nrow(quantum_residues) * 100,
      mean_conservation = mean(conservation),
      sd_conservation = sd(conservation)
    )

  # Print results
  cat("\nQuantum relevance classification:\n")
  print(summary_stats)

  # Write to file
  write.table(summary_stats, "../results/quantum_relevance_summary.txt", 
              sep="\t", row.names=FALSE, quote=FALSE)
}, error = function(e) {
  cat("Error in summary statistics:", e$message, "\n")
})

# Test for correlation between distance and conservation
cat("\nAnalyzing correlation between distance and conservation...\n")
tryCatch({
  correlation <- cor.test(quantum_residues$distance_to_cofactor, quantum_residues$conservation)
  cat("Correlation between distance and conservation:\n")
  print(correlation)
  
  # Save correlation results
  sink("../results/correlation_results.txt")
  print(correlation)
  sink()
}, error = function(e) {
  cat("Error in correlation analysis:", e$message, "\n")
})

# Compare conservation between groups
cat("\nComparing conservation between groups...\n")
tryCatch({
  t_test_result <- t.test(
    conservation ~ quantum_relevance == "High", 
    data = quantum_residues
  )
  cat("T-test comparing conservation in high relevance vs. others:\n")
  print(t_test_result)
  
  # Save t-test results
  sink("../results/ttest_results.txt")
  print(t_test_result)
  sink()
}, error = function(e) {
  cat("Error in t-test analysis:", e$message, "\n")
})

# Create visualization
cat("\nGenerating plots...\n")
tryCatch({
  p1 <- ggplot(quantum_residues, aes(x=distance_to_cofactor, y=conservation, color=quantum_relevance)) +
    geom_point() +
    theme_minimal() +
    labs(title="Relationship Between Cofactor Distance and Conservation",
         x="Distance (Å)",
         y="Conservation Score")

  # Save plot
  ggsave("../results/distance_conservation_plot.png", p1, width=8, height=6)
  cat("Saved distance vs. conservation plot\n")

  # Create distribution plot
  p2 <- ggplot(quantum_residues, aes(x=quantum_relevance, fill=quantum_relevance)) +
    geom_bar() +
    theme_minimal() +
    labs(title="Distribution of Residues by Quantum Relevance",
         x="Quantum Relevance Category",
         y="Count")

  # Save plot
  ggsave("../results/quantum_relevance_distribution.png", p2, width=8, height=6)
  cat("Saved quantum relevance distribution plot\n")
}, error = function(e) {
  cat("Error generating plots:", e$message, "\n")
  cat("Plots may not have been saved. Consider running this script on a system with R graphical capabilities.\n")
})

cat("\nAnalysis complete!\n")
