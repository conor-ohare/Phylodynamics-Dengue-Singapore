# Load libraries for phylogenetic analysis and data manipulation
library(tfpscanner)
library(ape)
library(tidytree)
library(fs)
library(dplyr)
library(lubridate)
library(purrr)
library(tools)

options(error = traceback)

# Retrieve command line arguments for project organization (expects 'time' and 'date')
args <- commandArgs(trailingOnly = TRUE)

# Set project directory based on the first command line argument
projectDir <- args[1]
setwd(paste0("/Users/conor/NUS/dengue/phylodynamics_dengue_singapore/", projectDir))

# Current date for naming conventions and file management
currentDate <- Sys.Date()

# Read phylogenetic metadata
metadata <- read.csv("tfpscanner_traits.csv")

# Initialize data frames
tfpGrowthRates <- data.frame(cluster_id = integer(0))
tfpCladeSizes <- data.frame(cluster_id = integer(0))

# List phylogenetic tree files in the 'trees' directory
treeFiles <- list.files(path = "./beast_trees/summary_trees")

# Ensure the 'plots' directory exists
if (!dir.exists("./beast_trees/tree_plots")) {
  dir.create("./beast_trees/tree_plots")}

# Loop through each tree file
for (treeFile in treeFiles) {
  # Read phylogenetic tree
  phyloTree <- read.nexus(paste0("./beast_trees/summary_trees/", treeFile))

  # Plot tree and save to tree_plots
  outputFilePath <- paste0("./beast_trees/tree_plots/", tools::file_path_sans_ext(treeFile), ".png")

  # Open a PNG device
  png(filename = outputFilePath)

  
  # Plot the tree to the PNG device
  plot(phyloTree, cex = 0.6, adj = 1)  # Adjust `cex` for font size and `adj` for text alignment

  
  # Close the PNG device
  dev.off()

  sampleIds <- phyloTree$tip.label
  relevantMetadata <- metadata[metadata$`sequence_name` %in% sampleIds, ]


  # Retrieve the sequence_name for the smallest sample_time
  root_name <- relevantMetadata$sequence_name[which.min(relevantMetadata$sample_time)]
  root_time <- relevantMetadata$sample_time[which.min(relevantMetadata$sample_time)]

  # Run tfp analysis with minimum descendants and reporting frequency
  tfpResults <- tfpscanner::tfpscan(tre = phyloTree, amd = relevantMetadata,
                                    min_descendants = 30, report_freq = 1,
                                    min_cluster_age_yrs = 0.1,
                                    root_on_tip = root_name,
                                    root_on_tip_sample_time = root_time,
                                    Tg = 90 / 365,
                                    detailed_output = TRUE, compute_gam = TRUE)



  # Generate a column name from the tree file name (sans extension)
  columnName <- file_path_sans_ext(basename(treeFile))

  # Prepare data frame with cluster IDs and calculated logistic growth rates
  tempDf_rate <- data.frame(cluster_id = tfpResults$cluster_id,
                       logistic_growth_rate = tfpResults$logistic_growth_rate)
  names(tempDf_rate) <- c("cluster_id", columnName)

  # Merge current results with the main tfp statistics data frame
  tfpGrowthRates <- merge(tfpGrowthRates, tempDf_rate, by = "cluster_id", all = TRUE)


  # Prepare data frame with cluster IDs and calculated logistic growth rates
  tempDf_size <- data.frame(cluster_id = tfpResults$cluster_id,
                            cluster_size = tfpResults$cluster_size)
  names(tempDf_size) <- c("cluster_id", columnName)

  # Merge current results with the main tfp statistics data frame
  tfpCladeSizes <- merge(tfpCladeSizes, tempDf_size, by = "cluster_id", all = TRUE)

}

# Save tfp growth rates statistics to CSV
write.csv(tfpGrowthRates, "tfp_growth_rates.csv", row.names = FALSE)

# Save tfp clade sizes statistics to CSV
write.csv(tfpCladeSizes, "tfp_clade_sizes.csv", row.names = FALSE)

## ----------------------------------------------------------------------------------------------------

## BROWSER INFORMATION
## important visualisation, but right now we're concerned with just generating
## the changing growth rates per cluster

## ----------------------------------------------------------------------------------------------------

# Define a path for the tfpscanner results, incorporating the current date for organization
# tfpscannerResultsDirName <- paste0("tfpscan-", currentDate)


# Move the tfpscanner results to the designated output directory
# file_move(tfpscannerResultsDirName, tfpscannerOutputDir)

# Create a directory for the tfpbrowser files within the tfpscanner output directory
  # tfpbrowserFilesDir <- file.path(tfpscannerOutputDir, "tfpbrowser_files")
  # dir_create(tfpbrowserFilesDir)

# # Prepare data for the tfpbrowser, specifying the environment file and output directory
# envFilePath <- file.path(tfpscannerOutputDir, tfpscannerResultsDirName, paste0("scanner-env-", currentDate, ".rds"))
# tfpscanner::create_browser_data(
#   e0 = envFilePath,
#   output_dir = tfpbrowserFilesDir
# )


# Open the generated HTML file in the default web browser to view the phylogenetic tree analysis
# htmlFilePath <- file.path(tfpbrowserFilesDir, "treeview/tree-logistic_growth_rate.html")
# browseURL(htmlFilePath)