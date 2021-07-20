
library(dplyr)
library(tidyverse)
library(janitor)
library(ggplot2)
library(patchwork)
library(airPLS)


setwd("C:/Users/kayte/Dropbox (MIT)/Tonegawa Lab/Projects/BLA Dopamine/Photometry/photometry_data_processing")
source("code/functions/movavg.R")
source("code/functions/get_zdFF.R")


RunPipeline <- function(files, path, FUN){
  out <- list()
  for(j in files){
    name <- str_remove(j, pattern = "^\\d+-")
    out[[paste(name)]] <- FUN(paste(path,j,sep = "/"))
  }
  return(out)
}

########### Clean Raw Data ##########################################################

# filename <- choose.files()
# df <- read_csv(filename, skip = 1)


input_folder <- "photometry_raw_data/photometry_cfc"
exp_name = "photometry"
file_list <- list.files(input_folder)


pipeline <- function(path) {
  df <- read_csv(filename, skip = 1)
  names(df) <- c("time", "ref_410nm", "da_470nm", "raw", "ttl1", "analog1", "analog2", "NA")
  start <- min(which(df$ttl1 == "1"))
  end <- max(which(df$ttl1 == "1"))
  df <- df[start:end,1:5]
  df <- mutate(df, time = (time - time[1]))
  raw_reference <- df$ref_410nm[1:nrow(df)]
  raw_signal <- df$da_470nm[1:nrow(df)]
  zdFF = get_zdFF(raw_reference, raw_signal, smooth_win= 60, lambda = 7e5)
  return(zdFF)
}

photometry_cfc <- RunPipeline(file_list, input_folder, pipeline)

# change column names
names(df) <- c("time", "ref_410nm", "da_470nm", "raw", "ttl1", "analog1", "analog2", "NA")

# Crop data to behavioral window, remove extra columns
start <- min(which(df$ttl1 == "1"))
end <- max(which(df$ttl1 == "1"))
df <- df[start:end,1:5]

# Reset time column
df <- mutate(df, time = (time - time[1]))

# # Plot raw data
# raw_ref_plot <- ggplot(df, aes(x= time, y =ref_410nm)) +
#   geom_line(colour="purple", size = .7) +
#   geom_vline(xintercept = c(180, 250, 320), size = 2, alpha = .3, colour = "red") +
#   xlab("Time(s)") +
#   ylab("Raw Signal") +
#   theme_minimal()
# 
# raw_da_plot <- ggplot(df, aes(x= time, y =da_470nm)) +
#   geom_line(colour="blue", size = .7) +
#   geom_vline(xintercept = c(180, 250, 320), size = 2, alpha = .3, colour = "red") +
#   xlab("Time(s)") +
#   ylab("Raw Signal") +
#   theme_minimal()
# 
# raw_ref_plot + raw_da_plot +
#   plot_layout(ncol=1, nrow = 2)

# Extract raw signals
raw_reference <- df$ref_410nm[1:nrow(df)]
raw_signal <- df$da_470nm[1:nrow(df)]

#### Run Analysis ####
zdFF = get_zdFF(raw_reference, raw_signal, smooth_win= 60, lambda = 7e5)
plot(zdFF, type='l', col='black')