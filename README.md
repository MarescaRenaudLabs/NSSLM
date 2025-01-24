**Data and code available at**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13988116.svg)](https://doi.org/10.5281/zenodo.13988116)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)

# NSSLM

This repository contains example code and data accompanying the paper: "Nonlinear sound-sheet microscopy: imaging opaque organs at the capillary and cellular scale" related to Nonlinear Sound Sheet Localization Microscopy (NSSLM).

## Folder Structure
- data: Folder containing the raw NSSLM data files downloaded from the Zenodo repository.
- display: Folder containing display-related functions
- examples: Folder containing example scripts for processing NSSLM data
  - 3DFFTFiltering: Subfolder containing an example script and images demonstrating how the 3D FFT filter works on a simple moving target
  - NSSLM: Subfolder containing an example script demonstrating how to perform NSSLM of the available dataset
- filtering: Folder containing the 3D FFT filtering function and SVD filtering function (from the PALA toolbox)

## Usage
### Example script: processNSSLM.m
This script demonstrates the post-processing of NSSM data to yield Non Linear Sound Sheet Localization Microscopy results.

**Key steps:**
1. Define Path List: Add necessary paths to the MATLAB environment.
2. Load Sequence Parameters: Load the sequence parameters from a .mat file.
3. Set Processing Parameters: Define various processing parameters such as filter velocity boundaries, SVD cutoff values, and more.
4. Load and Process Data: Load NSSM data files, apply pre-correlation, SVD filtering, and 3D FFT filtering.
5. Localization and Tracking: Perform localization and tracking of sound sheets using the PALA toolbox functions.
6. Save Results: Save the processed results and parameters to the specified directory.
7. Rendering: Calculate the size of each pixel, super-resolution image size, and adjust trajectories for rendering.

For more information see [readme file](examples/NSSLM/readme.md).

### Example script: 3DFFTfiltering.m
This script demonstrates space-time FFT filtering on a 2D target moving in the z or x-direction with different speeds

**Key steps**
1. Define Parameters: Set the size, pixel dimensions, PSF FWHM, target speed, and frame rate.
2. Create Target: Initialize a 3D matrix and create a meshgrid for the target.
3. Make Gaussian PSF: Generate a Gaussian Point Spread Function.
4. Simulate Movement: Update the target's position over time based on speed.
5. Apply FFT Filtering: Perform space-time FFT filtering on the generated sequence.

For more information see [readme file](examples/3DFFTfiltering/readme.md).

## Dependencies
- MATLAB
- PALA toolbox (available at [PALA GitHub](https://github.com/AChavignon/PALA))

## License
This project is licensed under the [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) License. See the LICENSE file for details.


## Acknowledgements
This project was created and maintained by Baptiste Heiles. Special thanks to the contributors and the PALA toolbox developers.
