library(tidyverse)
library(fs) # Great for file path handling

# --- 1. CONFIGURATION & FILE DISCOVERY ---

# Identify all CSV files in the current folder that look like results
# Pattern matches: ends in .csv
file_list <- list.files(pattern = "\\.csv$", full.names = FALSE)

# Filter out output files so we don't read them back in (infinite loop prevention)
file_list <- file_list[!str_detect(file_list, "Summary|calculated")]

if(length(file_list) == 0) stop("No input CSV files found in the working directory.")

cat(sprintf("Found %d files to process.\n", length(file_list)))

# --- 2. BATCH PROCESSING FUNCTION ---

process_single_file <- function(filepath) {
  
  # A. Read Data
  # check.names = FALSE ensures ImageJ column names like "Area" stay "Area"
  raw_data <- read_csv(filepath, show_col_types = FALSE)
  
  # B. Extract Metadata from Filename (Regex)
  f_name <- basename(filepath)
  
  # Extract PSI (allows decimals)
  psi_val <- str_extract(f_name, "([0-9.]+)(?=psi)") %>% as.numeric()
  
  # Extract Actuations (allows integers)
  act_val <- str_extract(f_name, "([0-9]+)(?=act)") %>% as.numeric()
  
  # Handle cases where regex fails (files named incorrectly)
  if(is.na(psi_val)) psi_val <- NA
  if(is.na(act_val)) act_val <- NA
  
  # C. Process Data
  processed_data <- raw_data %>%
    # Ensure we are using the correct column for diameter
    # ImageJ usually outputs 'Feret' for max diameter. 
    rename(diameter_mm = Feret) %>% 
    
    mutate(
      pressure_psi = psi_val,
      actuation_count = act_val,
      source_file = f_name,
      
      # CALCULATION: Sphere Volume = (Pi * d^3) / 6
      # Input: mm -> Output: mm^3 (which equals microliters, uL)
      volume_ul = (pi * (diameter_mm ^ 3)) / 6,
      
      # CONVERSION: uL -> nL (easier to read for your scale)
      volume_nl = volume_ul * 1000
    )
  
  return(processed_data)
}

# --- 3. EXECUTE BATCH LOOP ---

# map_dfr runs the function on every file and "binds" rows into one big table
master_dataset <- map_dfr(file_list, process_single_file)

# --- 4. CLEANUP & FORMATTING ---

final_table <- master_dataset %>%
  # Organize columns: Metadata first, then metrics
  select(
    pressure_psi,
    actuation_count,
    source_file,
    diameter_mm,
    volume_nl,  # Using nanoliters is likely better for your scale
    volume_ul
  ) %>%
  # Sort by actuation count so the table is readable
  arrange(pressure_psi, actuation_count) %>%
  # Round for display
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  # Add IDs
  mutate(Global_ID = row_number())

# --- 5. OUTPUT ---

# Print summary to console
print(head(final_table))

# Save the Master Dataset
write_csv(final_table, "Master_Droplet_Summary.csv")

cat("\nProcessing Complete. Saved to 'Master_Droplet_Summary.csv'.\n")



library(tidyverse)

# --- 1. CONFIGURATION ---
input_file <- "Master_Droplet_Summary.csv"

# Check if file exists before proceeding
if(!file.exists(input_file)) stop("Master_Droplet_Summary.csv not found. Run the previous batch script first.")

# --- 2. DATA LOADING ---
raw_data <- read_csv(input_file, show_col_types = FALSE)

# --- 3. STATISTICAL SUMMARY ---
# We calculate Mean, Standard Deviation (SD), and Coefficient of Variation (CV)
calibration_table <- raw_data %>%
  group_by(actuation_count) %>%
  summarise(
    avg_volume_nl = mean(volume_nl, na.rm = TRUE),
    sd_volume_nl  = sd(volume_nl, na.rm = TRUE),
    cv_percent    = (sd(volume_nl) / mean(volume_nl)) * 100,
    n_samples     = n()
  )

# Print the table to console so you can see the raw numbers immediately
print(calibration_table)

# Save the summary table
write_csv(calibration_table, "Calibration_Summary_Stats.csv")

# --- 4. LINEAR REGRESSION MODEL ---
# Ideally, Volume = Slope * Actuations + Intercept
model <- lm(avg_volume_nl ~ actuation_count, data = calibration_table)
model_summary <- summary(model)

# Extract stats for the plot label
r_sq <- round(model_summary$r.squared, 4)
slope <- round(coef(model)[2], 2)
intercept <- round(coef(model)[1], 2)

equation_label <- paste0("y = ", slope, "x + ", intercept, "\nR² = ", r_sq)

# --- 5. PLOTTING THE CURVE ---
cal_plot <- ggplot(calibration_table, aes(x = actuation_count, y = avg_volume_nl)) +
  
  # Error Bars (Standard Deviation) - This visualizes your chaotic compliance
  geom_errorbar(aes(ymin = avg_volume_nl - sd_volume_nl, 
                    ymax = avg_volume_nl + sd_volume_nl), 
                width = 0.2, color = "red") +
  
  # The Data Points
  geom_point(size = 4, color = "black") +
  
  # The Linear Fit Line
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  
  # Labels and Aesthetics
  labs(
    title = "Droplet Volume Calibration Curve",
    subtitle = paste("Linearity Check:", equation_label),
    x = "Number of Actuations (Pulse Count)",
    y = "Average Droplet Volume (nL)",
    caption = "Error bars represent Standard Deviation (SD)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

# Display the plot
print(cal_plot)

# Save the plot
ggsave("Calibration_Curve_Plot.png", plot = cal_plot, width = 6, height = 4, dpi = 300)

cat("\nAnalysis Complete. Check 'Calibration_Summary_Stats.csv' and 'Calibration_Curve_Plot.png'.\n")