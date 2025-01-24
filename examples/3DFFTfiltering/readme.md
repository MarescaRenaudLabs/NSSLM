# 3DFFTfiltering.m Usage Guide

## Overview
This script demonstrates space-time FFT filtering on a 2D target moving in the z-direction. It creates a synthetic target, applies a Gaussian Point Spread Function (PSF), and simulates the target's movement over time. The script then performs space-time FFT filtering to analyze the target's motion.

## Key Steps

### 1. Define Parameters
The script begins by defining the parameters for the target and the filtering process:
- `fkParam.size`: Size of the 3D target matrix [rows, columns, frames].
- `fkParam.sizeOfPixm`: Size of each pixel in meters [pixel_x, pixel_z].
- `fkParam.fwhmX`: Full width at half maximum (FWHM) of the Gaussian PSF in the x direction.
- `fkParam.fwhmZ`: Full width at half maximum (FWHM) of the Gaussian PSF in the z direction.
- `fkParam.speed`: Speed of the target in m/s [speed_x, speed_z].
- `fkParam.framerate`: Frame rate in Hz.

### 2. Create Target
The script initializes a 3D matrix `Target` with zeros and creates a meshgrid for the target. This meshgrid is used to define the spatial coordinates of the target.

### 3. Make Gaussian PSF
A Gaussian Point Spread Function (PSF) is created using the specified FWHM values. The PSF is used to simulate the target's appearance in the image.

### 4. Simulate Movement
The script simulates the movement of the target over time by updating its position based on the specified speed and frame rate. The PSF is applied to each frame to create a sequence of images representing the moving target.

### 5. Perform Space-Time FFT Filtering
The script performs space-time FFT filtering on the generated sequence of images to analyze the target's motion. This step involves applying a 3D FFT to the image sequence and filtering based on the target's velocity signature.

## Example 1: Target Moving in Z-Direction

### Parameters
- Size of the 3D target matrix: `[157, 161, 1000]`
- Size of each pixel in meters: `[56.6 56.6] * 1e-6`
- Full width at half maximum (FWHM) of the Gaussian PSF in x direction: `3`
- Full width at half maximum (FWHM) of the Gaussian PSF in z direction: `3`
- Speed of the target in m/s: `[15e-3, 0]`
- Frame rate in Hz: `1000`

### Create Target
- Initialize the target matrix with zeros.
- Create a meshgrid for the target.

### Make Gaussian PSF
- Generate a Gaussian Point Spread Function using the specified FWHM values.

### Simulate Movement
- Update the target's position over time based on the specified speed and frame rate.
- Apply the PSF to each frame to create a sequence of images representing the moving target.

### Perform Space-Time FFT Filtering
- Apply 3D FFT to the target matrix.
- Define the filter in the frequency domain.
- Apply the filter to the FFT of the target matrix.
- Perform inverse FFT to get the filtered target matrix.

### Filter and Display Positive Speeds
- Filter and display the filtered imaging sequence to retain only the positive speeds in the target range [5 20] mm/s.

### Filter and Display Negative Speeds
- Filter and display the filtered imaging sequence to retain only the negative speeds in the target range [5 20] mm/s.

## Example Images
### Original Target
![Original Target](figures/example1/target.gif)

### Filtered Target
#### Target filtered for negative speeds
![Filtered Target Negative Z Speeds](figures/example1/targetFilteredNeg.gif)
#### Target filtered for positive speeds
![Filtered Target Positive Z Speeds](figures/example1/targetFilteredPos.gif)

## Example 2: Target Moving in X-Direction

### Define Parameters
- Size of the 3D target matrix (e.g. [157, 161, 1000])  
- Size of each pixel in meters (e.g. [56.6 56.6] * 1e-6)  
- FWHM (e.g. 3 in both x and z directions)  
- Speed of the target in m/s (e.g. [0, 15e-3])  
- Frame rate in Hz (e.g. 1000)

### Create Target
- Initialize the target matrix with zeros.  
- Create a meshgrid for the target.  
- Apply the Gaussian PSF at each frame.

### Perform Space-Time FFT Filtering
- Apply 3D FFT to the target matrix.  
- Define the filter in the frequency domain (e.g. based on speed range).  
- Apply the filter and perform inverse FFT to get the filtered target matrix.

### Filter and Display Positive Speeds
- Filter to retain only positive speeds in a specified range.  
- Display the resulting imaging sequence.

### Filter and Display Negative Speeds
- Filter to retain only negative speeds in that range.  
- Display the resulting imaging sequence.

## Example Images

### Original Target
![Original Target](figures/example2/target.gif)

### Filtered Target
#### Target Filtered for Negative X Speeds
![Filtered Target Negative X Speeds](figures/example2/targetFilteredNeg.gif)  
#### Target Filtered for Positive X Speeds
![Filtered Target Positive X Speeds](figures/example2/targetFilteredPos.gif)

---

## Example 3: Bidirectional Motion (Z and X)

### Define Parameters
- Size of the 3D target matrix (e.g. [157, 161, 1000])  
- Size of each pixel in meters (e.g. [56.6 56.6] * 1e-6)  
- FWHM (e.g. 3 in both x and z directions)  
- Speed of the target in m/s (e.g. [10e-3, 10e-3])  
- Frame rate in Hz (e.g. 1000)

### Create Target
- Initialize the target matrix with zeros.  
- Create a meshgrid for the target.  
- Apply the Gaussian PSF at each frame, updating both Z and X positions.

### Perform Space-Time FFT Filtering
- Apply 3D FFT to the sequence.  
- Filter based on your chosen velocity bounds for both Z and X directions.  
- Inverse FFT to obtain the filtered sequence.

### Filter and display the filtered imaging sequence to retain only the x-positive speeds in target range [5 20] mm/s
- Focus on positive velocity signatures in the selected range for X only.  
- Display the sequence.

### Filter and display the filtered imaging sequence to retain the z-negative and x-positive speeds in target range [5 20] mm/s
- Focus on positive velocity signatures in the selected range for Z and X only.  
- Display the sequence.

### Filter and display the filtered imaging sequence to retain the x and z-negative speeds in target range [5 20] mm/s
- Focus on negative velocity signatures in the selected range for Z and all speeds for both and X.  
- Display the sequence.

### Filter and display the filtered imaging sequence to retain only the x-negative speeds in target range [5 20] mm/s
- Focus on positive velocity signatures in the selected range for Z and all speeds for both and X.  
- Display the sequence.

## Example Images

### Original Target
![Original Target](figures/example3/target.gif)

### Filtered Target
#### Target Filtered for x-positive speeds in target range [5 20] mm/s
![Target Filtered for x-positive speeds in target range [5 20] mm/s](figures/example3/targetFilteredPosx.gif)

#### Target Filtered for z-negative and x-positive speeds in target range [5 20] mm/s
![Target Filtered for z-negative and x-positive speeds in target range [5 20] mm/s](figures/example3/targetFilteredNegzPosx.gif)

#### Target Filtered for x and z-negative speeds in target range [5 20] mm/s
![Target Filtered for x and z-negative speeds in target range [5 20] mm/s](figures/example3/targetFilteredNegzNegx.gif)

#### Target Filtered for x-negative speeds in target range [5 20] mm/s
![Target Filtered for x-negative speeds in target range [5 20] mm/s](figures/example3/targetFilteredNegx.gif)


## Example 4: Multiple targets and speeds

### Define Parameters
- Size of the 3D target matrix (e.g. [157, 161, 1200])  
- Pixel size in meters (e.g. [56.6 56.6] * 1e-6)  
- FWHM (e.g. 3 in both x and z directions)  
- A more complex speed profile (e.g. [0, 15e-3] changing over time, or multiple speed phases)  
- Frame rate in Hz (e.g. 1000)

### Create Targets
- Initialize the target matrix with zeros.  
- Create a meshgrid.  
- Simulate a motion profile that changes speed over time, or multiple objects moving at different speeds.

### Perform Space-Time FFT Filtering
- Apply 3D FFT to the entire sequence.  
- Design filters to isolate specific velocity ranges or complex velocity patterns.  
- Inverse FFT to reconstruct the filtered sequence.

### Filter and Display x positive speeds in target range [5-20] mm/s
- Apply the filter to capture x-positive velocities in the chosen range(s).  
- Display the filtered sequence.

### Filter and Display x positive speeds in target range [15-20] mm/s
- Apply the filter to capture x-positive velocities in the chosen range(s).  
- Compare with the positive-velocity sequence to analyze the bidirectional or multi-speed motion.

## Example Images

### Original Target
![Original Target](figures/example4/target.gif)

### Filtered Target
#### Target Filtered for  x positive speeds in target range [5-20] mm/s
![Target Filtered for  x positive speeds in target range [5-20] mm/s](figures/example4/targetFilteredPosx_5_20.gif)

#### Target Filtered for  x positive speeds in target range [15-20] mm/s
![Target Filtered for  x positive speeds in target range [15-20] mm/s](figures/example4/targetFilteredPosx_15_20.gif)