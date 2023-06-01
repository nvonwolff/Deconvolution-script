import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.cm as cm

#WORK ONGOING
# Load the observed spectrum and reference spectra from Excel files
# ---Make sure that the files are in Excel format and that the first column holds the wavelength and the second the absorption data---
observed_spectrum_file = 'observed_spectrum.xlsx'
# ---Add the name of as many reference spectra as you want. Make sure that the script is in the same folder as the reference and observed spectra---
reference_spectra_files = ['reference_spectrum1.xlsx', 'reference_spectrum2.xlsx', 'reference_spectrum2.xlsx']

# Define the wavelength range (in nm) for truncation
wavelength_range = (230, 600)

# Load observed spectrum and truncate observed spectrum
observed_spectrum_df = pd.read_excel(observed_spectrum_file)
observed_wavelength = observed_spectrum_df.iloc[:, 0].values
observed_intensity = observed_spectrum_df.iloc[:, 1].values
observed_mask = np.logical_and(observed_wavelength >= wavelength_range[0], observed_wavelength <= wavelength_range[1])
observed_wavelength = observed_wavelength[observed_mask]
observed_intensity = observed_intensity[observed_mask]

# Initialize arrays to store truncated and interpolated reference spectra
reference_wavelengths = []
reference_intensities = []

# Load and truncate reference spectra
for file in reference_spectra_files:
    reference_spectrum_df = pd.read_excel(file)
    reference_wavelength = reference_spectrum_df.iloc[:, 0].values
    reference_intensity = reference_spectrum_df.iloc[:, 1].values

    # Truncate reference spectrum
    reference_mask = np.logical_and(reference_wavelength >= wavelength_range[0], reference_wavelength <= wavelength_range[1])
    reference_wavelength = reference_wavelength[reference_mask]
    reference_intensity = reference_intensity[reference_mask]

    reference_wavelengths.append(reference_wavelength)
    reference_intensities.append(reference_intensity)

# Find the minimum number of data points among all spectra within the defined wavelength range
min_num_points = min(len(observed_wavelength), *[len(reference_wavelength) for reference_wavelength in reference_wavelengths])

# Randomly delete excess data points from the observed spectrum
if len(observed_wavelength) > min_num_points:
    delete_indices = np.random.choice(len(observed_wavelength), len(observed_wavelength) - min_num_points, replace=False)
    observed_wavelength = np.delete(observed_wavelength, delete_indices)
    observed_intensity = np.delete(observed_intensity, delete_indices)

# Randomly delete excess data points from the reference spectra
for i in range(len(reference_wavelengths)):
    num_points = len(reference_wavelengths[i])
    if num_points > min_num_points:
        delete_indices = np.random.choice(num_points, num_points - min_num_points, replace=False)
        reference_wavelengths[i] = np.delete(reference_wavelengths[i], delete_indices)
        reference_intensities[i] = np.delete(reference_intensities[i], delete_indices)

# Normalize all spectra (Unmark if you want to activate this feature)
# observed_intensity /= np.max(observed_intensity)
# for i in range(len(reference_intensities)):
#     reference_intensities[i] /= np.max(reference_intensities[i])

# Define the fitting function as a linear combination of reference spectra
def fit_func(wavelength, *contributions):
    fitted_spectrum = np.zeros_like(wavelength)
    for i, contribution in enumerate(contributions):
        fitted_spectrum += contribution * reference_intensities[i]
    return fitted_spectrum

# Initial guess for the contributions
initial_contributions = np.ones(len(reference_intensities))

# Perform curve fitting to estimate the contributions
popt, _ = curve_fit(fit_func, observed_wavelength, observed_intensity, p0=initial_contributions)

# Calculate the fitted spectrum using the estimated contributions
fitted_spectrum = fit_func(observed_wavelength, *popt)

# Calculate the residuals
residuals = observed_intensity - fitted_spectrum

# Calculate R-squared value
ss_total = np.sum((observed_intensity - np.mean(observed_intensity)) ** 2)
ss_residual = np.sum(residuals ** 2)
r_squared = 1 - (ss_residual / ss_total)

# Calculate the percentages of the reference spectra contributions
total_intensity = np.sum(popt)
percentages = [contribution / total_intensity * 100 for contribution in popt]

# Determine the number of reference spectra
num_reference_spectra = len(reference_intensities)

# Generate a colormap with different shades of blue
cmap = cm.get_cmap('YlGn')
colors = cmap(np.linspace(0.4, 0.8, num_reference_spectra))

# Plot the observed spectrum, fitted spectrum, and contributions
plt.figure(figsize=(6, 6))
plt.scatter(observed_wavelength, observed_intensity, s=30, facecolors='none', edgecolors='r',
            label=f'Observed {observed_spectrum_file.split(".")[0]}')
plt.plot(observed_wavelength, fitted_spectrum, color='black', alpha=0.75, linewidth=4,
         label=f'Fitted Spectrum, R2 = {r_squared:.4f}')

for i, contribution in enumerate(popt):
    color = colors[i]
    plt.plot(observed_wavelength, contribution * reference_intensities[i], '--', color=color, linewidth=3,
             label=f'weighted {reference_spectra_files[i].split(".")[0]} ({percentages[i]:.2f}%)')

plt.xlabel('Wavelength (nm)', fontweight='bold')
plt.ylabel('Absorption (a.u.)', fontweight='bold')
plt.xticks(fontweight='bold')
plt.yticks(fontweight='bold')
plt.legend(frameon=False)

# Save the figure as a PNG file
output_filename = observed_spectrum_file.split(".")[0] + ".png"
plt.savefig(output_filename, dpi=300)
print(f"Figure saved as {output_filename}")

# Show the plot
plt.tight_layout()
plt.show()

# Print the percentages
for i, percentage in enumerate(percentages):
    print(f'Percentage of Contribution {i + 1}: {percentage:.2f}%')
    
# Save the fitted spectrum to an Excel file
def save_fitted_spectrum_to_excel(output_file):
    fitted_spectrum_data = pd.DataFrame({'Wavelength': observed_wavelength, 'Intensity': fitted_spectrum})
    fitted_spectrum_data.to_excel(output_file, index=False)
    print(f"Fitted spectrum saved to {output_file}")

# Specify the output filename for the fitted spectrum (change name if necessary)
output_spectrum_file = 'FittedSpectrum.xlsx'

# Call the function to save the fitted spectrum
save_fitted_spectrum_to_excel(output_spectrum_file)
