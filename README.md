# Spectra-Deconvolution
Small code to deconvolute and fit UVVis spectra created with the help of chatGPT

# Spectral Fitting Script

This script fits an observed spectrum with a linear combination of reference spectra. It utilizes the Python libraries `pandas`, `numpy`, `matplotlib`, and `scipy` for data manipulation, visualization, and curve fitting. It thus allows to calculate the amount of different species contributing to the spectroscopical response of a mixture. This can be used for example to track the built-up or consumption of species over time. 

## Prerequisites

- Python 3.x
- pandas
- numpy
- matplotlib
- scipy

## Branch usage

Changes should be made to the fixes branch.

## Usage

1. Ensure that the script file is located in the same directory as the observed and reference spectra files, which should be in Excel format.
2. Update the script with the appropriate file names for the observed spectrum and reference spectra. Make sure that the first column of each file contains the wavelength values, and the second column contains the absorption data.
3. Optionally, modify the wavelength range to truncate the spectra according to your needs.
4. Optionally, uncomment the normalization code if you want to normalize the spectra.
5. Run the script.

The script will perform the following steps:

1. Load and truncate the observed spectrum and reference spectra within the specified wavelength range.
2. If necessary, randomly delete excess data points from the observed spectrum and reference spectra to ensure equal lengths.
3. Fit the observed spectrum using a linear combination of the reference spectra.
4. Calculate the residuals and R-squared value to evaluate the fit quality.
5. Generate a plot showing the observed spectrum, fitted spectrum, and contributions from each reference spectrum.
6. Save the plot as a PNG file.
7. Print the percentages of the reference spectra contributions.
8. Saves the fitted spectrum as an excel file

## Additional Customizations

- You can modify the colormap used for the contributions plot by changing the `cmap` value in the script. Currently, it uses the 'YlGn' colormap.
- Adjust the color shades of the contributions plot by modifying the `colors` array. The script currently generates different shades of green.

Feel free to customize the script according to your specific requirements.

# Notes
The reference spectra and possibly the observed spectrum should be corrected for the concentration of species before running the script.

## License

This script is released under the [MIT License](https://opensource.org/licenses/MIT).
