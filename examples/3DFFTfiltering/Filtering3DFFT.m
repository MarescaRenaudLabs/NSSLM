%% This script demonstrates a simple example of space-time FFT filtering on a simple target moving in 2D
% Created by Baptiste Heiles on 2025/01/13
% Last modified by Baptiste Heiles on 2025/01/22

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% First example with a target moving in z-direction %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
fkParam.size = [157, 161, 1000]; % Size of the 3D target matrix [rows, columns, frames]
fkParam.sizeOfPixm = [56.6 56.6] * 1e-6; % Size of each pixel in meters [pixelZ, pixelX]
fkParam.fwhmX = 3; % Full width at half maximum (FWHM) of the Gaussian PSF in x direction
fkParam.fwhmZ = 3; % Full width at half maximum (FWHM) of the Gaussian PSF in z direction
fkParam.speed = [15e-3, 0]; % Speed of the target in m/s [speedZ, speedX]
fkParam.framerate = 1000; % Frame rate in Hz

%% Create target
Target = zeros(fkParam.size); % Initialize the target matrix with zeros
szz = floor(fkParam.size(1) / 2); % Half size in z direction
szx = floor(fkParam.size(2) / 2); % Half size in x direction
[meshTgtX, meshTgtZ] = meshgrid(-szx:szx, -szz:szz); % Create meshgrid for target
meshTgtIn = cat(3, meshTgtX, meshTgtZ); % Concatenate meshgrid into a 3D array

%% Make Gaussian PSF
vectfwhmX = -fkParam.fwhmX:fkParam.fwhmX; % Vector for FWHM in x direction
vectfwhmZ = -fkParam.fwhmZ:fkParam.fwhmZ; % Vector for FWHM in z direction
[meshX, meshZ] = meshgrid(vectfwhmX, vectfwhmZ); % Create meshgrid for PSF
meshIn = cat(3, meshX, meshZ); % Concatenate meshgrid into a 3D array
sigGaussZ = vectfwhmZ(end) * 0 + 1; % Standard deviation for Gaussian in z direction
sigGaussX = vectfwhmX(end) * 0 + 1; % Standard deviation for Gaussian in x direction
myGaussFunc = @(posX, meshPos) (exp(- (meshPos(:, :, 1) - posX(1)) .^ 2 / (2 * sigGaussZ ^ 2) - (meshPos(:, :, 2) - posX(2)) .^ 2 / (2 * sigGaussX ^ 2))); % Gaussian function
posX = [0 0]; % Initial position for PSF
PSF = exp(- (meshIn(:, :, 1) - posX(1)) .^ 2 / (2 * sigGaussZ ^ 2) - (meshIn(:, :, 2) - posX(2)) .^ 2 / (2 * sigGaussX ^ 2)); % Create Gaussian PSF

%% Create all frames with translating PSF based on speed
posXInit = [-szx + 45, -szz + 20]; % Initial position of the target [x, z]
posXUpd = posXInit; % Updated position of the target
stepZ = (fkParam.speed(1) / fkParam.framerate) / fkParam.sizeOfPixm(1); % Step size in z direction per frame
stepX = (fkParam.speed(2) / fkParam.framerate) / fkParam.sizeOfPixm(2); % Step size in x direction per frame

for iFr = 1:fkParam.size(3)
    posXUpd = posXUpd + [stepX, stepZ]; % Update position of the target

    % Handle boundary conditions for x direction
    if posXUpd(1) > szx
        posXUpd(1) = -szx;
        stepX = -stepX;
    elseif posXUpd(1) < -szx
        posXUpd(1) = szx;
        stepX = -stepX;
    end

    % Handle boundary conditions for z direction
    if posXUpd(2) > szz
        posXUpd(2) = szz;
        stepZ = -stepZ;
    elseif posXUpd(2) < -szz
        posXUpd(2) = -szz;
        stepZ = -stepZ;
    end

    % Add the Gaussian PSF to the target matrix at the updated position
    Target(:, :, iFr) = Target(:, :, iFr) + myGaussFunc(posXUpd, meshTgtIn);
end

% Display the target matrix using a logarithmic scale
displayImg(Target, 'log', 'none');

%% Filter and display the filtered imaging sequence to retain only the positive speeds in target range [5 20] mm/s
maxSpeed = vecnorm(fkParam.speed, 2); % Calculate the maximum speed of the target
[TargetFilt, fftImgFilt] = filt3DFFT(Target, [5e-3 20e-3], [fkParam.sizeOfPixm 1 ./ fkParam.framerate], 1, 0);
displayImg(TargetFilt, 'log', 'none', [-25 0]);

%% Filter and display the filtered imaging sequence to retain only the negative speeds in target range [5 20] mm/s
maxSpeed = vecnorm(fkParam.speed, 2); % Calculate the maximum speed of the target
[TargetFilt, fftImgFilt] = filt3DFFT(Target, [5e-3 20e-3], [fkParam.sizeOfPixm 1 ./ fkParam.framerate], -1, 0);
displayImg(TargetFilt, 'log', 'none', [-25 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Second example with a target moving in x-direction %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
fkParam.size = [157, 161, 1000]; % Size of the 3D target matrix [rows, columns, frames]
fkParam.sizeOfPixm = [56.6 56.6] * 1e-6; % Size of each pixel in meters [pixelZ, pixelX]
fkParam.fwhmX = 3; % Full width at half maximum (FWHM) of the Gaussian PSF in x direction
fkParam.fwhmZ = 3; % Full width at half maximum (FWHM) of the Gaussian PSF in z direction
fkParam.speed = [0, 15e-3]; % Speed of the target in m/s [speedZ, speedX]
fkParam.framerate = 1000; % Frame rate in Hz

%% Create target
Target = zeros(fkParam.size); % Initialize the target matrix with zeros
szz = floor(fkParam.size(1) / 2); % Half size in z direction
szx = floor(fkParam.size(2) / 2); % Half size in x direction
[meshTgtX, meshTgtZ] = meshgrid(-szx:szx, -szz:szz); % Create meshgrid for target
meshTgtIn = cat(3, meshTgtX, meshTgtZ); % Concatenate meshgrid into a 3D array

%% Make Gaussian PSF
vectfwhmX = -fkParam.fwhmX:fkParam.fwhmX; % Vector for FWHM in x direction
vectfwhmZ = -fkParam.fwhmZ:fkParam.fwhmZ; % Vector for FWHM in z direction
[meshX, meshZ] = meshgrid(vectfwhmX, vectfwhmZ); % Create meshgrid for PSF
meshIn = cat(3, meshX, meshZ); % Concatenate meshgrid into a 3D array
sigGaussZ = vectfwhmZ(end) * 0 + 1; % Standard deviation for Gaussian in z direction
sigGaussX = vectfwhmX(end) * 0 + 1; % Standard deviation for Gaussian in x direction
myGaussFunc = @(posX, meshPos) (exp(- (meshPos(:, :, 1) - posX(1)) .^ 2 / (2 * sigGaussZ ^ 2) - (meshPos(:, :, 2) - posX(2)) .^ 2 / (2 * sigGaussX ^ 2))); % Gaussian function
posX = [0 0]; % Initial position for PSF
PSF = exp(- (meshIn(:, :, 1) - posX(1)) .^ 2 / (2 * sigGaussZ ^ 2) - (meshIn(:, :, 2) - posX(2)) .^ 2 / (2 * sigGaussX ^ 2)); % Create Gaussian PSF

%% Create all frames with translating PSF based on speed
posXInit = [-szx + 45, -szz + 20]; % Initial position of the target [x, z]
posXUpd = posXInit; % Updated position of the target
stepZ = (fkParam.speed(1) / fkParam.framerate) / fkParam.sizeOfPixm(1); % Step size in z direction per frame
stepX = (fkParam.speed(2) / fkParam.framerate) / fkParam.sizeOfPixm(2); % Step size in x direction per frame

for iFr = 1:fkParam.size(3)
    posXUpd = posXUpd + [stepX, stepZ]; % Update position of the target

    % Handle boundary conditions for x direction
    if posXUpd(1) > szx
        posXUpd(1) = szx;
        stepX = -stepX;
    elseif posXUpd(1) < -szx
        posXUpd(1) = -szx;
        stepX = -stepX;
    end

    % Handle boundary conditions for z direction
    if posXUpd(2) > szz
        posXUpd(2) = szz;
        stepZ = -stepZ;
    elseif posXUpd(2) < -szz
        posXUpd(2) = -szz;
        stepZ = -stepZ;
    end

    % Add the Gaussian PSF to the target matrix at the updated position
    Target(:, :, iFr) = Target(:, :, iFr) + myGaussFunc(posXUpd, meshTgtIn);
end

% Display the target matrix using a logarithmic scale
displayImg(Target, 'log', 'none');

%% Filter and display the filtered imaging sequence to retain only the x-positive speeds in target range [5 20] mm/s
maxSpeed = vecnorm(fkParam.speed, 2); % Calculate the maximum speed of the target
[TargetFilt, fftImgFilt] = filt3DFFT(Target, [5e-3 20e-3], [fkParam.sizeOfPixm 1 ./ fkParam.framerate], 0, 1);
displayImg(TargetFilt, 'log', 'none', [-25 0]);

%% Filter and display the filtered imaging sequence to retain only the x-negative speeds in target range [5 20] mm/s
maxSpeed = vecnorm(fkParam.speed, 2); % Calculate the maximum speed of the target
[TargetFilt, fftImgFilt] = filt3DFFT(Target, [5e-3 20e-3], [fkParam.sizeOfPixm 1 ./ fkParam.framerate], 0, -1);
displayImg(TargetFilt, 'log', 'none', [-25 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Third example with a target moving in z-x-direction %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
fkParam.size = [157, 161, 1000]; % Size of the 3D target matrix [rows, columns, frames]
fkParam.sizeOfPixm = [56.6 56.6] * 1e-6; % Size of each pixel in meters [pixelZ, pixelX]
fkParam.fwhmX = 3; % Full width at half maximum (FWHM) of the Gaussian PSF in x direction
fkParam.fwhmZ = 3; % Full width at half maximum (FWHM) of the Gaussian PSF in z direction
fkParam.speed = 15e-3 .* [1, 1]; % Speed of the target in m/s [speedZ, speedX]
fkParam.framerate = 1000; % Frame rate in Hz

%% Create target
Target = zeros(fkParam.size); % Initialize the target matrix with zeros
szz = floor(fkParam.size(1) / 2); % Half size in z direction
szx = floor(fkParam.size(2) / 2); % Half size in x direction
[meshTgtX, meshTgtZ] = meshgrid(-szx:szx, -szz:szz); % Create meshgrid for target
meshTgtIn = cat(3, meshTgtX, meshTgtZ); % Concatenate meshgrid into a 3D array

%% Make Gaussian PSF
vectfwhmX = -fkParam.fwhmX:fkParam.fwhmX; % Vector for FWHM in x direction
vectfwhmZ = -fkParam.fwhmZ:fkParam.fwhmZ; % Vector for FWHM in z direction
[meshX, meshZ] = meshgrid(vectfwhmX, vectfwhmZ); % Create meshgrid for PSF
meshIn = cat(3, meshX, meshZ); % Concatenate meshgrid into a 3D array
sigGaussZ = vectfwhmZ(end) * 0 + 1; % Standard deviation for Gaussian in z direction
sigGaussX = vectfwhmX(end) * 0 + 1; % Standard deviation for Gaussian in x direction
myGaussFunc = @(posX, meshPos) (exp(- (meshPos(:, :, 1) - posX(1)) .^ 2 / (2 * sigGaussZ ^ 2) - (meshPos(:, :, 2) - posX(2)) .^ 2 / (2 * sigGaussX ^ 2))); % Gaussian function
posX = [0 0]; % Initial position for PSF
PSF = exp(- (meshIn(:, :, 1) - posX(1)) .^ 2 / (2 * sigGaussZ ^ 2) - (meshIn(:, :, 2) - posX(2)) .^ 2 / (2 * sigGaussX ^ 2)); % Create Gaussian PSF

%% Create all frames with translating PSF based on speed
posXInit = [0, -szz]; % Initial position of the target [x, z]
posXUpd = posXInit; % Updated position of the target
stepZ = (fkParam.speed(1) / fkParam.framerate) / fkParam.sizeOfPixm(1); % Step size in z direction per frame
stepX = (fkParam.speed(2) / fkParam.framerate) / fkParam.sizeOfPixm(2); % Step size in x direction per frame

for iFr = 1:fkParam.size(3)
    posXUpd = posXUpd + [stepX, stepZ]; % Update position of the target

    % Handle boundary conditions for x direction
    if posXUpd(1) > szx
        posXUpd(1) = szx;
        stepX = -stepX;
    elseif posXUpd(1) < -szx
        posXUpd(1) = -szx;
        stepX = -stepX;
    end

    % Handle boundary conditions for z direction
    if posXUpd(2) > szz
        posXUpd(2) = szz;
        stepZ = -stepZ;
    elseif posXUpd(2) < -szz
        posXUpd(2) = -szz;
        stepZ = -stepZ;
    end

    % Add the Gaussian PSF to the target matrix at the updated position
    Target(:, :, iFr) = Target(:, :, iFr) + myGaussFunc(posXUpd, meshTgtIn);
end

% Display the target matrix using a logarithmic scale
displayImg(Target, 'log', 'none');

%% Filter and display the filtered imaging sequence to retain only the x-positive speeds in target range [5 20] mm/s
maxSpeed = vecnorm(fkParam.speed, 2); % Calculate the maximum speed of the target
[TargetFilt, fftImgFilt] = filt3DFFT(Target, [5e-3 20e-3], [fkParam.sizeOfPixm 1 ./ fkParam.framerate], 0, 1);
displayImg(TargetFilt, 'log', 'none', [-25 0]);

%% Filter and display the filtered imaging sequence to retain the z-negative and x-positive speeds in target range [5 20] mm/s
maxSpeed = vecnorm(fkParam.speed, 2); % Calculate the maximum speed of the target
[TargetFilt, fftImgFilt] = filt3DFFT(Target, [5e-3 20e-3], [fkParam.sizeOfPixm 1 ./ fkParam.framerate], -1, 1);
displayImg(TargetFilt, 'log', 'none', [-25 0]);

%% Filter and display the filtered imaging sequence to retain the x and z-negative speeds in target range [5 20] mm/s
maxSpeed = vecnorm(fkParam.speed, 2); % Calculate the maximum speed of the target
[TargetFilt, fftImgFilt] = filt3DFFT(Target, [5e-3 20e-3], [fkParam.sizeOfPixm 1 ./ fkParam.framerate], -1, -1);
displayImg(TargetFilt, 'log', 'none', [-25 0]);

%% Filter and display the filtered imaging sequence to retain only the negative speeds in target range [5 20] mm/s
maxSpeed = vecnorm(fkParam.speed, 2); % Calculate the maximum speed of the target
[TargetFilt, fftImgFilt] = filt3DFFT(Target, [5e-3 20e-3], [fkParam.sizeOfPixm 1 ./ fkParam.framerate], 0, -1);
displayImg(TargetFilt, 'log', 'none', [-25 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Fourth example with multiple targets and speeds %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
fkParam.size = [157, 161, 1000]; % Size of the 3D target matrix [rows, columns, frames]
fkParam.sizeOfPixm = [56.6 56.6] * 1e-6; % Size of each pixel in meters [pixelZ, pixelX]
fkParam.fwhmX = 3; % Full width at half maximum (FWHM) of the Gaussian PSF in x direction
fkParam.fwhmZ = 3; % Full width at half maximum (FWHM) of the Gaussian PSF in z direction
fkParam.speed = 50e-3 .* [1, 1]; % Speed of the target in m/s [speedZ, speedX]
fkParam.framerate = 1000; % Frame rate in Hz

%% Create target
Target = zeros(fkParam.size); % Initialize the target matrix with zeros
szz = floor(fkParam.size(1) / 2); % Half size in z direction
szx = floor(fkParam.size(2) / 2); % Half size in x direction
[meshTgtX, meshTgtZ] = meshgrid(-szx:szx, -szz:szz); % Create meshgrid for target
meshTgtIn = cat(3, meshTgtX, meshTgtZ); % Concatenate meshgrid into a 3D array

%% Make Gaussian PSF
vectfwhmX = -fkParam.fwhmX:fkParam.fwhmX; % Vector for FWHM in x direction
vectfwhmZ = -fkParam.fwhmZ:fkParam.fwhmZ; % Vector for FWHM in z direction
[meshX, meshZ] = meshgrid(vectfwhmX, vectfwhmZ); % Create meshgrid for PSF
meshIn = cat(3, meshX, meshZ); % Concatenate meshgrid into a 3D array
sigGaussZ = vectfwhmZ(end) * 0 + 1; % Standard deviation for Gaussian in z direction
sigGaussX = vectfwhmX(end) * 0 + 1; % Standard deviation for Gaussian in x direction
myGaussFunc = @(posX, meshPos) (exp(- (meshPos(:, :, 1) - posX(1)) .^ 2 / (2 * sigGaussZ ^ 2) - (meshPos(:, :, 2) - posX(2)) .^ 2 / (2 * sigGaussX ^ 2))); % Gaussian function

%% Create several targets
for idx_target = 1:4
    posX = [0 0]; % Initial position for PSF
    PSF = exp(- (meshIn(:, :, 1) - posX(1)) .^ 2 / (2 * sigGaussZ ^ 2) - (meshIn(:, :, 2) - posX(2)) .^ 2 / (2 * sigGaussX ^ 2)); % Create Gaussian PSF

    %% Create all frames with translating PSF based on speed
    posXInit = idx_target .* 10 .* [1, 0] + [0, -szz]; % Initial position of the target [x, z]
    posXUpd = posXInit; % Updated position of the target
    stepZ = 1 ./ idx_target .* (fkParam.speed(1) / fkParam.framerate) / fkParam.sizeOfPixm(1); % Step size in z direction per frame
    stepX = 1 ./ idx_target .* (fkParam.speed(2) / fkParam.framerate) / fkParam.sizeOfPixm(2); % Step size in x direction per frame

    for iFr = 1:fkParam.size(3)
        posXUpd = posXUpd + [stepX, stepZ]; % Update position of the target

        % Handle boundary conditions for x direction
        if posXUpd(1) > szx
            posXUpd(1) = szx;
            stepX = -stepX;
        elseif posXUpd(1) < -szx
            posXUpd(1) = -szx;
            stepX = -stepX;
        end

        % Handle boundary conditions for z direction
        if posXUpd(2) > szz
            posXUpd(2) = szz;
            stepZ = -stepZ;
        elseif posXUpd(2) < -szz
            posXUpd(2) = -szz;
            stepZ = -stepZ;
        end

        % Add the Gaussian PSF to the target matrix at the updated position
        Target(:, :, iFr) = Target(:, :, iFr) + myGaussFunc(posXUpd, meshTgtIn);
    end

end

% Display the target matrix using a logarithmic scale
displayImg(Target, 'log', 'none');

%% Filter and display the filtered imaging sequence to retain only the x-positive speeds in target range [5 20] mm/s
maxSpeed = vecnorm(fkParam.speed, 2); % Calculate the maximum speed of the target
[TargetFilt, fftImgFilt] = filt3DFFT(Target, [5e-3 20e-3], [fkParam.sizeOfPixm 1 ./ fkParam.framerate], 0, 1);
displayImg(TargetFilt, 'log', 'none', [-25 0]);

%% Filter and display the filtered imaging sequence to retain the x-positive speeds in target range [15 20] mm/s
maxSpeed = vecnorm(fkParam.speed, 2); % Calculate the maximum speed of the target
[TargetFilt, fftImgFilt] = filt3DFFT(Target, [15e-3 20e-3], [fkParam.sizeOfPixm 1 ./ fkParam.framerate], 0, 1);
displayImg(TargetFilt, 'log', 'none', [-25 0]);
