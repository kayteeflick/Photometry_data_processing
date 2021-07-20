
library(dplyr)
library(tidyverse)
library(janitor)
library(ggplot2)
library(patchwork)

library(airPLS)
source("code/functions/movavg.R")
source("code/functions/get_zdFF.R")

########### Clean Raw Data ##########################################################

filename <- choose.files()
df <- read_csv(filename, skip = 1)

# change column names
names(df) <- c("time", "ref_410nm", "da_470nm", "raw", "ttl1", "analog1", "analog2", "NA")

# Crop data to behavioral window, remove extra columns
start <- min(which(df$ttl1 == "1"))
end <- max(which(df$ttl1 == "1"))
df <- df[start:end,1:5]

# Reset time column
df <- mutate(df, time = (time - time[1]))

# Plot raw data
raw_ref_plot <- ggplot(df, aes(x= time, y =ref_410nm)) +
  geom_line(colour="purple", size = .7) +
  geom_vline(xintercept = c(180, 250, 320), size = 2, alpha = .3, colour = "red") +
  xlab("Time(s)") +
  ylab("Raw Signal") +
  theme_minimal()

raw_da_plot <- ggplot(df, aes(x= time, y =da_470nm)) +
  geom_line(colour="blue", size = .7) +
  geom_vline(xintercept = c(180, 250, 320), size = 2, alpha = .3, colour = "red") +
  xlab("Time(s)") +
  ylab("Raw Signal") +
  theme_minimal()

raw_ref_plot + raw_da_plot +
  plot_layout(ncol=1, nrow = 2)

# Extract raw signals
raw_reference <- df$ref_410nm[1:nrow(df)]
raw_signal <- df$da_470nm[1:nrow(df)]

########## Baseline Subtraction ###############################################

### Smooth

# Run function movavg.R
smooth_win = 24 #24 corresponds to 200ms
smooth_reference <- movavg(raw_reference, smooth_win)
smooth_signal <- movavg(raw_signal, smooth_win)

# Plot smoothed signal
par(mfrow=c(3,1))
plot(smooth_reference, type='l', col='purple')
plot(smooth_signal, type='l', col='blue')

### Find slope baseline

# Install package airPLS 
# install.packages('devtools')
# library(devtools)
# httr::set_config( httr::config( ssl_verifypeer = 0L ) )
# install_github("zmzhang/airPLS_R")

# remove bleaching curve with iterative matching

#inv_ref <-  1/smooth_reference
#base_r <- airPLS(inv_ref,5e5,1,20)
#base_r <- 1/base_r
base_r <- airPLS(smooth_reference,5e5,1,20)
base_r_4 <- airPLS(smooth_reference,5e4,1,20)
base_r_6 <- airPLS(smooth_reference,5e6,1,20)
base_s <- airPLS(smooth_signal,5e5,1,20)
base_s_4 <- airPLS(smooth_signal,5e4,1,20)
base_s_6 <- airPLS(smooth_signal,5e6,1,20)

# Plot data and baselines
par(mfrow=c(2,1))
ref_baseline_plot <- plot(x = df$time, y = smooth_reference, 
                          type='l', col='purple', xlab = "Time(s)", ylab = "Reference")
lines(x = df$time, y = base_r, col = "orange", )
lines(x = df$time, y = base_r_4, col = "black")
lines(x = df$time, y = base_r_6, col = "red")
da_baseline_plot <- plot(x = df$time,y = smooth_signal, 
                         type='l', col='blue', xlab = "Time(s)", ylab = "Signal")
lines(x = df$time, y = base_s, col = "orange")
lines(x = df$time, y = base_s_4, col = "black")
lines(x = df$time, y = base_s_6, col = "red")


### Remove baseline

reference <- raw_reference - base_r
signal <- raw_signal - base_s

# Plot flatten data
par(mfrow=c(2,1))
plot(reference, type='l', col='purple')
plot(signal, type='l', col='blue')

##################### PART 2 #####################

### Calculate dF/F

movavg_win =200 #120 corresponds to 1s
movavg_reference <- movavg(reference, movavg_win)
movavg_signal <- movavg(signal, movavg_win)

plot(movavg_reference, type='l', col='purple')
plot(movavg_signal, type='l', col='blue')

reference_dff <-  (reference - movavg_reference)/movavg_reference
signal_dff <- (signal - movavg_signal)/movavg_signal

# Plot dF/F data
par(mfrow=c(2,1))
plot(reference_dff, type='l', col='purple')
plot(signal_dff, type='l', col='blue')


### Standardize signals
z_reference <- (reference - median(reference)) / sd(reference)
z_signal <- (signal - median(signal)) / sd(signal)

# Plot standardized data
par(mfrow=c(2,1))
plot(z_reference, type='l', col='purple')
plot(z_signal, type='l', col='blue')

### Linear robust fit
require(MASS)
fit <- rlm(z_signal ~ z_reference)

# Plot fit
par(mfrow=c(1,1))
plot(reference, signal) 
abline(fit, col="red")


### Align signals

z_reference_fit <- predict(fit)

# Plot aligned signals
plot(z_signal, type='l', col='blue')
lines(z_reference_fit, type='l', col='purple')


### Calculate z-score dF/F

# subtract reference channel
zdFF <- z_signal - z_reference_fit

# Plot z-score dF/F
plot(zdFF, type='l', col='black')

