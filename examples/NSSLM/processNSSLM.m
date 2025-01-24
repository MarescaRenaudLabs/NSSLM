%% A script for post-processing of NSSM data to yield Non Linear Sound Sheet Localization Microscopy
% Created by Baptiste Heiles on 2024/12/25
% Last modified by Baptiste Heiles on 2025/01/22

%% Define path list
Paths.pathToMain = fullfile('..', '..'); % Add path to main directory containing the different repos NSSLM and PALA
Paths.pathToRepo = fullfile(Paths.pathToMain); % Define path to the NSSLM repository
Paths.pathToPALA = fullfile(Paths.pathToMain, 'PALA', 'PALA', 'PALA_addons'); % Define path to PALA toolbox addons (for more information, see https://github.com/AChavignon/PALA)
Paths.pathToData = fullfile(Paths.pathToMain, 'data'); % Define path to the data folder

addpath(fullfile(Paths.pathToRepo, 'filtering')); % Add all subfolders in the NSSLM repository to the path
addpath(fullfile(Paths.pathToRepo, 'display'));
addpath(genpath(Paths.pathToPALA)); % Add all subfolders in the PALA addons to the path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% NSSLM processing pipeline %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define directories and saving prefixes
Paths.folderToSave = fullfile(Paths.pathToData, 'Soundsheet_1'); % Define folder to save results
Paths.filePrefix = 'ImgNSSM'; % Define file prefix for NSSM files
numFiles = size(dir(fullfile(Paths.pathToData, [Paths.filePrefix '*.mat'])), 1); % Get number of NSSM files
Paths.filename = sprintf('%s_%03i.mat', Paths.filePrefix, 1); % Define filename for the first NSSM file
Paths.folderToFFT = fullfile(Paths.folderToSave, 'fft_Img_test'); % Define folder to save FFT images
Paths.filePrefixFFT = 'fft_Img'; % Define file prefix for FFT images

%% Example with a filter in between 3-15 mm/s for Soundsheet 1
PPParam.iSS = 1; % Set index for Soundsheet. Data is 4D with sound-sheet in 3rd dimension

PPParam.is3DFFTfilt = 1; % Enable 3D FFT filtering
PPParam.filtVelocity = [0 3] .* 1e-3; % Set filter velocity boundaries in m/s
PPParam.isPreCorr = 1; % Enable pre-filter 1-lag correlation
PPParam.isSVDfilt = 1; % Enable SVD filtering
PPParam.SVDCutoff = [2, 2100]; % Set SVD cutoff values
PPParam.filterSign = [0 0]; % sign of filter on kz kx
timestamp = datestr(now, 'yyyy_mm_dd_HH_MM_ss'); % Get current timestamp
Paths.pathToResult = fullfile(Paths.folderToSave, sprintf('Result_%1.1f_%1.1f_mms_%s', PPParam.filtVelocity(1) .* 1e3, PPParam.filtVelocity(2) .* 1e3, timestamp)); % Define path to save results
mkdir(Paths.pathToResult); % Create directory to save results

%% Compute scaling factors, wavelengths etc
PPParam.TransmitFrequency = 13600000; % Transmit frequency in Hz
PPParam.speedOfSound = 1540; % Speed of sound in m/s
PPParam.lambda = PPParam.speedOfSound ./ PPParam.TransmitFrequency .* 1e3; % Wavelength in mm
PPParam.framerate = 1000; % Frame rate in Hz
PPParam.units = [1/2 1/2 1 / PPParam.framerate] .* [PPParam.lambda .* 1e-3 PPParam.lambda .* 1e-3 1]; % Define units

%% Define rendering parameters
IntPower = 1/5;
SigmaGauss = 3;

%% Load pre-filtered NSSM images
tProc = tic; % Start timer

for iFiles = 1:numFiles
    Paths.fileName = sprintf('%s_%03d.mat', Paths.filePrefix, iFiles); % Define filename for current NSSM file
    % Load mat
    Img = importdata(fullfile(Paths.pathToData, Paths.fileName));

    % Crop Img
    croppedImg = squeeze(Img(:, :, PPParam.iSS, :));

    % 1-lag correlation for signal enhancement
    if PPParam.isPreCorr
        ImgCorr = croppedImg(:, :, [1:end - 1]) .* conj(croppedImg(:, :, [2:end]));
    else
        ImgCorr = croppedImg;
    end

    % SVD filter
    if PPParam.isSVDfilt
        ImgFiltSVD = SVDfilter(ImgCorr, PPParam.SVDCutoff);
    end

    % 3D FFT filter
    if PPParam.is3DFFTfilt
        [Img3DFFTFilt, FFTImgFilt, ~] = filt3DFFT(ImgFiltSVD, PPParam.filtVelocity, PPParam.units, PPParam.filterSign(1), PPParam.filterSign(2));
        % For convenience, you can save fftImgFilt once instead of
        % computing it everytime for each of the velocity bands you want

        % Save fft_Img_filt
        % mkdir(fullfile(Paths.folderToFFT));
        % save(fullfile(Paths.folderToFFT, sprintf('%s_%03d.mat', Paths.filePrefixFFT, iFiles)), 'FFTImgFilt');

    else
        Img3DFFTFilt = ImgFiltSVD;
    end

    % The next steps use the ULM functions available in the PALA toolbox
    % (available at: https://github.com/AChavignon/PALA)

    % Define the ULM structure
    res = 10;
    ULM = struct('numberOfParticles', 80, ... % Number of particles per frame. (30-100)
        'size', [size(Img3DFFTFilt)], ... % size of input data [nb_pixel_z nb_pixel_x nb_frame_per_bloc]
        'scale', [1/2 1/2 1 / PPParam.framerate], ... % Scale [z x dt], size of pixel in the scaling unit. (here, pixsize = 1*lambda)
        'res', res, ... % Resolution factor. Typically 10 for final image rendering at lambda/10.
        'SVD_cutoff', [1 size(Img3DFFTFilt)], ... % SVD filtering, to be adapted to your clutter/SNR levels
        'max_linking_distance', 1, ... % Maximum linking distance between two frames to reject pairing, in pixels units (UF.scale(1)). (2-4 pixel).
        'min_length', 3, ... % Minimum allowed length of the tracks in time. (5-20 frames)
        'fwhm', [1 1] * 3, ... % Size [pixel] of the mask for localization. (3x3 for pixel at lambda, 5x5 at lambda/2). [fmwhz fmwhx]
        'max_gap_closing', 0, ... % Allowed gap in microbubbles' pairing. (if you want to skip frames 0)
        'interp_factor', 1 / res, ... % Interpfactor (decimation of tracks)
        'LocMethod', 'Radial' ... % Select localization algorithm (WA,Interp,Radial,CurveFitting,NoLocalization)
    );
    ULM.parameters.NLocalMax = 3; % Safeguard on the number of maxLocal in the fwhm*fwhm grid (3 for fwhm=3, 7 for fwhm=5)
    ULM.lambda = PPParam.lambda; scale = ULM.scale;

    % Detection and localization process (return a list of coordinates in pixel)
    [Localizations{iFiles}] = ULM_localization2D(abs(Img3DFFTFilt), ULM);

    % Tracking algorithm (list of tracks)
    Track_tot_i = ULM_tracking2D(Localizations{iFiles}, ULM);
    Trajectories{iFiles} = Track_tot_i;

end

tend = toc(tProc); disp('Done')
fprintf('ULM done in %d hours %.1f minutes.\n', floor(tend / 60/60), rem(tend / 60, 60));
mkdir(Paths.pathToResult);
tSave = tic;

% Save the data in different cells
mkdir(fullfile(Paths.pathToResult, 'Localizations'));
mkdir(fullfile(Paths.pathToResult, 'Trajectories'));
clear LocalizationsCell TrajectoriesCell

for iCell = 1:size(Localizations, 2)
    LocalizationsCell = Localizations{iCell};
    TrajectoriesCell = Trajectories{iCell};

    save(fullfile(Paths.pathToResult, 'Localizations', sprintf('LocalizationsCell_%03i.mat', iCell)), 'LocalizationsCell');
    save(fullfile(Paths.pathToResult, 'Trajectories', sprintf('TrajectoriesCell_%03i.mat', iCell)), 'TrajectoriesCell');
end

save(fullfile(Paths.pathToResult, 'Params.mat'), 'ULM', 'PPParam');
tendSave = toc(tSave);
fprintf('Saving done in %d hours %.1f minutes.\n', floor(tendSave / 60/60), rem(tendSave / 60, 60));

%% Example of a rendering (similarly to what is described in PALA toolbox)
TrajectoriesVec = cat(1, Trajectories{:}); % Put all cells for one file in column
ULM.SizeOfPixelZm = ULM.lambda .* ULM.scale(1) .* 1e-3; % Define size of pixel in meters
ULM.SizeOfPixelXm = ULM.lambda .* ULM.scale(2) .* 1e-3; % Define size of pixel in meters
ULM.SizeOfPixelZmRendering = 1e-6; %pixe size in meters
ULM.SizeOfPixelXmRendering = 1e-6; %pixe size in meters
ULM.SRsize = round(ULM.size(1:2) .* [ULM.SizeOfPixelZm ULM.SizeOfPixelXm] ./ [ULM.SizeOfPixelZmRendering ULM.SizeOfPixelZmRendering]);
ULM.res_bypass = [ULM.SizeOfPixelZm ./ ULM.SizeOfPixelZmRendering, ULM.SizeOfPixelXm ./ ULM.SizeOfPixelXmRendering];
TrajectoriesHybrid = cellfun(@(x)((x(:, [1 2 3 4 5]) - [1, 1, 1, 1, 0]) .* [ULM.res_bypass(1), ULM.res_bypass(2), ULM.SizeOfPixelZm, ULM.SizeOfPixelXm, 1]), TrajectoriesVec, 'UniformOutput', 0); % Store all trajectories into one single matrix
MatrixTrack = ULM_Track2MatOut(TrajectoriesHybrid, ULM.SRsize(1:2) + [1 1] * 1, 'mode', '2D_tracks'); % build  matrix

%% Final Render ULM image (for a more detailed description, see ULM toolbox and accompanyin renderingNSSLM.m file)
% Intensity power and Gaussian sigma for rendering
% IntPower and SigmaGauss are defined earlier in the script
ULM_Image = MatrixTrack .^ IntPower;
if SigmaGauss > 0; ULM_Image = imgaussfilt(ULM_Image,SigmaGauss);end
    
figure('Position', [600 150 900 800]);
zax = [0:ULM.SRsize(1)] * ULM.SizeOfPixelZmRendering;
xax = [0:ULM.SRsize(2)] * ULM.SizeOfPixelXmRendering;

imagesc(xax * 1e3, zax * 1e3, ULM_Image); % % Display image
colormap(hot); % Set color map to hot
axis image; % Set axis to image
title(sprintf('ULM intensity map velocities: %s mm/s', mat2str(PPParam.filtVelocity * 1e3, 2)))
colorbar
% Note: The generated image has less tracks then Fig 5E, since its is
%       based on 6300 frames instead of 100000.
