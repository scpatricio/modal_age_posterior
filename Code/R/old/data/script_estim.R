setwd("/Users/scp93/Desktop/PhD/codes/mode_dist/data/")

# Cleaning the work enviroment 
rm(list = ls())

# Close any open connections (optional, but good practice)
closeAllConnections()

# Load required packages
library(dplyr)
library(parallel)
library(tidyr)


# Source required functions
source("/Users/scp93/Desktop/PhD/codes/mode_dist/mode_exact.R")

countries = c("US", "France", "Italy", "Japan", "Denmark", "Netherlands")
age_0 = 40
country = countries[1]



mode_summary = data.frame(Country = NULL,
                          Sex = NULL,
                          Year = NULL,
                          mean = NULL,
                          mode = NULL,
                          sd_value = NULL,
                          min_mode = NULL,
                          max_mode = NULL)

for(country in countries){
  # loading death counts data
  file_name_dx = paste0(country, "_Deaths_1x1.txt")
  df_dx = read.table(file = file_name_dx, 
                     header = TRUE,            # Use header because columns are named
                     sep = "",                 # Use whitespace as separator
                     skip = 2,                 # Skip the first two lines
                     stringsAsFactors = FALSE) # Keep text columns as character type
  
  df_new_dx = df_dx %>% 
    mutate(Age = as.numeric(Age),
           Age = ifelse(is.na(Age), 110, Age)) %>% 
    filter(!is.na(Age)) %>%
    pivot_longer(cols = c(Female, Male, Total), 
                 names_to = "Sex", 
                 values_to = "Deaths") %>% 
    filter(Year <= 2021 & Year >= 1960, 
           Sex != "Total")

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

  # Display the first few rows of the merged dataframe
  head(df_merged)
  
  
  sexes = unique(df_new_dx$Sex)
  years = unique(df_new_dx$Year)
  for(sex in sexes){
    for(year in years){
      print(paste(country, sex, year))
      # Filter data for current sex and year
      df_aux <- df_merged %>% 
        filter(Sex == sex, 
               Year == year, 
               Age >= age_0,
               Exposures > 0)
      
      # Estimate mode distribution
      mode_estim <- mode_dist_gauss(x = df_aux$Age, dx = df_aux$Deaths, ex = df_aux$Exposures)
      x = df_aux$Age
      dx = df_aux$Deaths
      ex = df_aux$Exposures
      
      # Calculate the mean of the distribution
      mean_mode = sum(mode_estim$mode * mode_estim$prob, na.rm = T)
      # Calculate the median of the distribution
      mode_mode = mean(mode_estim$mode[which.max(mode_estim$prob)], na.rm = T)
      # Calculate the variance
      variance <- sum((mode_estim$mode - mean_mode)^2 * mode_estim$prob, na.rm = T)
      # Calculate the standard deviation
      sd_value <- sqrt(variance)
      
      idx = which(mode_estim$prob>0)
      
      if(length(idx)>0){
        # Calculate the minimum of the distribution
        min_mode = min(mode_estim$mode[idx], na.rm = T)
        # Calculate the maximum of the distribution
        max_mode = max(mode_estim$mode[idx], na.rm = T)
      } else {
        min_mode = max_mode = NA
      }
      
      # Append summary data
      dx_aux = data.frame(
        Country = country,
        Sex = sex,
        Year = year,
        mean = mean_mode,
        mode = mode_mode,
        sd_value = sd_value,
        min_mode = min_mode,
        max_mode = max_mode
      )
      
      mode_summary = rbind(mode_summary, dx_aux)
    }
  }
}