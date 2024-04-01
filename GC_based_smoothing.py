import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Load the data from the CSV file
df = pd.read_csv('MFE_1042353N_GAPDH201.csv')

# Drop rows containing NaNs or Infs in 'GC' or 'GC_sumif' columns
df = df.replace([np.inf, -np.inf], np.nan)  # Replace all Infs with NaNs
df = df.dropna(subset=['GC_category', 'GC_sumif'])  # Drop rows with NaNs in specified columns

GC_fragment = df['GC_category'].values
fragment_count = df['GC_sumif'].values

# The function model for fitting
def quad_func(x, a, b, c, d):
    return a + b * x + c * np.power(x, 2) + d * np.power(x, 3)

# Attempt curve fitting on the cleaned data
params, params_covariance = curve_fit(quad_func, GC_fragment, fragment_count)

# Definition for using fitted parameters
def smoothed_function(x):
    return quad_func(x, *params)

# Add predictions back into DataFrame (optional step depending on your goal)
df['predicted_GC_sumif'] = df['GC_category'].apply(smoothed_function)
df['predicted_GC_sumif'] = df['predicted_GC_sumif'].clip(lower=0)  # Ensure non-negative predictions

# Save to CSV (as described in previous instructions)
output_filename = 'smooth_data_with_predictions_GC.csv'
df.to_csv(output_filename, index=False)

print(f"Output saved to {output_filename}")
