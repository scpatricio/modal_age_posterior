# Set working directory and clean environment
setwd("/Users/scp93/Desktop/PhD/codes/mode_dist/data/")
rm(list = ls())

# Close any open connections (optional, but good practice)
closeAllConnections()

# Load required packages
library(dplyr)
library(tidyr)
library(parallel)

# Source required functions
source("/Users/scp93/Desktop/PhD/codes/mode_dist/mode_exact.R")

# Define parameters
countries <- c("US", "France", "Italy", "Japan", "Denmark", "Netherlands")
age_0 <- 40

# Pre-allocate lists to store results
mode_summary_list <- list()
mode_full_sample_list <- list()

# Parallel processing setup
cl <- makeCluster(detectCores() - 1)  # Use one less than the available cores
clusterExport(cl, c("mode_dist_gauss", "age_0"))  # Export needed variables to clusters
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(tidyr))
clusterEvalQ(cl, library(mvtnorm))

# Loop over each country in parallel
results <- parLapply(cl, countries, function(country) {
  
  # Read country data
  file_name <- paste0(country, "_Deaths_1x1.txt")
  df_dx <- read.table(file = file_name, header = TRUE, sep = "", skip = 2, stringsAsFactors = FALSE)
  
  # Data transformation and filtering
  df_new_dx <- df_dx %>% 
    mutate(Age = as.numeric(Age),
           Age = ifelse(is.na(Age), 110, Age)) %>% 
    pivot_longer(cols = c(Female, Male, Total), names_to = "Sex", values_to = "Deaths") %>% 
    filter(!is.na(Age), Year <= 2021 & Year >= 1960)
  
  # loading exposure data
  file_name_ex = paste0(country, "_Exposures_1x1.txt")
  df_ex = read.table(file = file_name_ex,
                     header = TRUE,            # Use header because columns are named
                     sep = "",                 # Use whitespace as separator
                     skip = 2,                 # Skip the first two lines
                     stringsAsFactors = FALSE) # Keep text columns as character type
  
  df_new_ex = df_ex %>%
    mutate(Age = as.numeric(Age),
           Age = ifelse(is.na(Age), 110, Age)) %>%
    filter(!is.na(Age)) %>%
    pivot_longer(cols = c(Female, Male, Total),
                 names_to = "Sex",
                 values_to = "Exposures") %>%
    filter(Year <= 2021 & Year >= 1960,
           Sex != "Total")
  
  # Merge df_new_ex and df_new_dx by Year, Sex, and Age
  df_merged <- left_join(df_new_dx, df_new_ex, by = c("Year", "Sex", "Age"), )
  
  # Initialize lists to store intermediate results
  country_summary <- list()
  country_full_sample <- list()
  
  # Loop over sex and year combinations
  sexes <- unique(df_merged$Sex)
  years <- unique(df_merged$Year)
  for(sex in sexes) {
    for(year in years) {
      # Filter data for current sex and year
      df_aux <- df_merged %>% 
        filter(Sex == sex, Year == year, Age >= age_0, Exposures * Deaths > 0)
      
      if(nrow(df_aux) > 0){
        # Estimate mode distribution
        mode_estim <- mode_dist_gauss(x = df_aux$Age, dx = df_aux$Deaths, ex = df_aux$Exposures)
        
        # Calculate the mean of the distribution
        mean_mode = sum(mode_estim$mode * mode_estim$prob, na.rm = T)
        # Calculate the median of the distribution
        mode_mode = mean(mode_estim$mode[which.max(mode_estim$prob)], na.rm = T)
        # Calculate the variance
        variance <- sum((mode_estim$mode - mean_mode)^2 * mode_estim$prob, na.rm = T)
        # Calculate the standard deviation
        sd_value <- sqrt(variance)
        
        # Append summary data
        country_summary[[length(country_summary) + 1]] <- data.frame(
          Country = country,
          Sex = sex,
          Year = year,
          mean = mean_mode,
          mode = mode_mode,
          sd_value = sd_value,
          quantile_2_5 = mode_estim$quantile_2_5[1],
          quantile_97_5 = mode_estim$quantile_97_5[1]
        )
      }
    }
  }
  
  # Combine lists into data frames
  list(
    summary = do.call(rbind, country_summary),
    full_sample = do.call(rbind, country_full_sample)
  )
})

# Stop parallel cluster
stopCluster(cl)

# Combine results from all countries
mode_summary <- do.call(rbind, lapply(results, `[[`, "summary"))
mode_full_sample <- do.call(rbind, lapply(results, `[[`, "full_sample"))

# Output the results

write.table(mode_summary, "mode_summary_new.txt", row.names = F)
write.table(mode_full_sample, "mode_full_sample_new.txt", row.names = F)
head(mode_summary)
# mode_full_sample
