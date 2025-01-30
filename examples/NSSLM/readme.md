# processNSSLM.m

## Overview
This MATLAB script demonstrates the post-processing of NSSM (Non-Linear Sound Sheet Microscopy) data to produce super-resolved images and trajectories of sound sheets. It incorporates pre-correlation, SVD filtering, and 3D FFT filtering, before performing a localization and tracking step.

## Requirements
- MATLAB
- PALA toolbox (https://github.com/AChavignon/PALA)  
- Required MATLAB functions in the "functions" folder

## Key Steps
1. **Define Paths**  
   Sets up all repository paths and loads required data, including sequence parameters.
2. **Data Loading**  
   Loads 4D NSSM data and performs cropping to isolate a specified sound sheet.
3. **Pre-Correlation**  
   Optionally applies a one-lag correlation to improve signal quality.
4. **SVD Filtering**  
   Applies SVD-based clutter filtering to remove noise and low-frequency components.
5. **3D FFT Filtering**  
   Filters the data based on a velocity range in the space-time domain.
6. **Localization & Tracking**  
   Uses ULM functions in the PALA toolbox to detect and track microbubble over time.
7. **Saving**  
   Results are saved to separate folders, including localizations, trajectories, and computed parameters.
8. **Rendering (Optional)**  
   A final section illustrates the super-resolution rendering steps.

## Usage
1. Clone or download this repository from Zenodo.
2. Update paths in the script (e.g., `Paths.pathToData`, `Paths.pathToRepo`) to match your local directory structure.  
3. Run the script in MATLAB:
   ```matlab
   processNSSLM
   ```
## Data 
Example data is available on Zenodo [zenodo.14740878](https://doi.org/zenodo.14740878). 

| File                                          | Content                                                                                                                                                                                                                     |
| --------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ImgNSSM_00[1-3].mat`                         | Reconstructed Nonlinear Sound Sheets Microscopy images. 3 datasets, each containing 2 sound sheets and 2100 frames per file. Data is part of the dataset to create Fig 5E. Process to a NSSLM image using `processNSSLM.m`. |
| `Trajectories/TrajectoriesCell_[001-050].mat` | Trajectories obtained after running `processNSSLM.m` on the full recording of Fig 5E. Load trajectories and use `renderingNSSLM.m` to visualize.                                                                            |
| `Params.mat`                                  | Acquisition and Processing parameters to create Fig 5E.                                                                                                                                                                     |