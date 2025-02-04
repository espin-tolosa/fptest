import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# Reading the CSV filename from argv
if len(sys.argv) < 2:
    print("Usage: python pearson_correlation.py <filename>")
    sys.exit(1)

filename = sys.argv[1]  # Read the first argument

# Function to format arrays as C-style
def format_c_array(name, array):
    formatted_values = ",\n".join(f"{val:.23e}" for val in array)
    return f"{formatted_values},"

def compute_coefficients(x):
    x_m = np.mean(x)
    sum_diff = np.sum(x - x_m)
    return x_m, sum_diff

def load_csv_data(filename):
    data = pd.read_csv(filename, header=None)  # Read CSV without headers
    x = data.iloc[:, 0].to_numpy()  # First column as NumPy array
    y = data.iloc[:, 1].to_numpy()  # Second column as NumPy array
    return x, y

# Compute Pearson correlation coefficient
# Create an array with 100 elements
x, y = load_csv_data( filename )

corr_coefficient, _ = pearsonr(x, y)

mean_val, sum_differences = compute_coefficients(x)

#print("P: ", mean_val, sum_differences)

# Print the arrays in C format

# Write the arrays to files
#with open("x_data.txt", "w") as fx:
#    fx.write(format_c_array("x", x))

#with open("y_data.txt", "w") as fy:
#    fy.write(format_c_array("y", y))

# Write Pearson correlation coefficient to a file
with open("pearson_correlation.txt", "w") as fc:
    fc.write(f"{corr_coefficient:.23e}\n")

#print(f"Expected: {corr_coefficient:.23e}\n")

#print(f"Sum of x: {np.sum(x):.16f}")
#print(f"Sum of y: {np.sum(y):.16f}")
#print(f"Sum of x*x: {np.sum(x*x):.16f}")
#print(f"Sum of y*y: {np.sum(y*y):.16f}")
#print(f"Sum of x*y: {np.sum(x*y):.16f}")
