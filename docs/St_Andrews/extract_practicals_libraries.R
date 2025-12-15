#!/usr/bin/env Rscript
# Extract libraries from practicals folder and save as R data file

extract_practicals_libraries <- function() {
  practicals_dir <- "practicals"
  
  if (!dir.exists(practicals_dir)) {
    stop("Practicals folder not found: ", practicals_dir)
  }
  
  # Get all .qmd files from practicals folder
  qmd_files <- list.files(
    path = practicals_dir, 
    pattern = "\\.qmd$", 
    full.names = TRUE,
    ignore.case = TRUE,
    recursive = TRUE  # Include subfolders if any
  )
  
  cat("Found", length(qmd_files), "QMD files in practicals folder\n")
  
  all_libraries <- character(0)
  lib_usage <- list()
  
  for (file in qmd_files) {
    cat("Scanning:", file, "\n")
    content <- tryCatch(
      readLines(file, warn = FALSE),
      error = function(e) character(0)
    )
    
    # Find library calls
    lib_lines <- grep("^library\\(", content, value = TRUE)
    if (length(lib_lines) > 0) {
      lib_names <- gsub("^library\\(([^)]+)\\).*", "\\1", lib_lines)
      lib_names <- gsub('["\']', '', lib_names)
      lib_names <- trimws(lib_names)
      
      all_libraries <- c(all_libraries, lib_names)
      
      # Track which files use which libraries
      for (lib in lib_names) {
        if (lib %in% names(lib_usage)) {
          lib_usage[[lib]] <- c(lib_usage[[lib]], basename(file))
        } else {
          lib_usage[[lib]] <- basename(file)
        }
      }
    }
  }
  
  # Get unique libraries
  unique_libs <- unique(sort(all_libraries))
  
  # Create summary data
  libraries_summary <- data.frame(
    Library = unique_libs,
    UsedIn = sapply(unique_libs, function(lib) {
      files <- unique(lib_usage[[lib]])
      paste(files, collapse = ", ")
    }),
    stringsAsFactors = FALSE
  )
  
  # Save the data
  saveRDS(libraries_summary, file = "practicals_libraries.rds")
  
  # Also save as CSV for easy viewing
  write.csv(libraries_summary, "practicals_libraries.csv", row.names = FALSE)
  
  cat("✅ Saved library data:\n")
  cat("   - practicals_libraries.rds (for Quarto)\n")
  cat("   - practicals_libraries.csv (for viewing)\n")
  cat("📊 Found", length(unique_libs), "unique libraries across", length(qmd_files), "files\n")
  
  return(libraries_summary)
}

# Run the function
extract_practicals_libraries()