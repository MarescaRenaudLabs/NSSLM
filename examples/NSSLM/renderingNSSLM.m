%% A script for rendering of NSSLM trajectory data
% Created by Baptiste Heiles on 2024/12/25
% Last modified by Baptiste Heiles on 2025/01/22

%% Define path list
Paths.pathToMain = fullfile('..', '..'); % Add path to main directory containing the different repos NSSLM and PALA
Paths.pathToRepo = fullfile(Paths.pathToMain); % Define path to the NSSLM repository
Paths.pathToData = fullfile(Paths.pathToMain, 'data'); % Define path to the data folder
Paths.pathToTrajectories = fullfile(Paths.pathToData, 'Trajectories'); % Define path to the directory holding all the TrajectoriesCell data .mat files

addpath(fullfile(Paths.pathToRepo, 'display'));

%% Setup Rendering
% Load the Params.mat matrix
load(fullfile(Paths.pathToData, 'Params.mat'));

% Define saving directory
saveFolder = fullfile(Paths.pathToData, 'Renderings'); % Define path to the directory where images will be saved
saveModeImg = 1; % Save images (1) or not (0)
timestamp = datestr(now, 'yyyy_mm_dd_HH_MM_ss');

% Define the scale of the units of the loaded trajectories and of the final image (for more information see PALA toolbox)
ULM.SizeOfPixelZm = ULM.lambda .* ULM.scale(1) .* 1e-3; % Scale of Z in meters
ULM.SizeOfPixelXm = ULM.lambda .* ULM.scale(2) .* 1e-3; % Scale of X in meters
ULM.SRsize = ULM.size(1:2) .* ULM.res; % SRsize is the size of the rendered image in pixels

%% Load all trajectories
TrajectoriesAll = [];
for is = 1:50
    tmpTrajectories = importdata(fullfile(Paths.pathToTrajectories, sprintf('TrajectoriesCell_%03i.mat', is)));
    if isempty(TrajectoriesAll)
        TrajectoriesAll = tmpTrajectories;
    else
        TrajectoriesAll = cat(1, TrajectoriesAll, tmpTrajectories);
    end
end
clear tmpTrajectories

% Conver the Trajectories to resolution factor and to meters
TrajectoriesMeters = [];
TrajectoriesResFact = cellfun(@(x)((x(:, [1 2 3 4 5]) - [1, 1, 0, 0, 0]) .* [ULM.res, ULM.res, 1, 1, 1]), TrajectoriesAll, 'UniformOutput', 0); % Store the trajectories with the resolution factor included
TrajectoriesMeters = cellfun(@(x)((x(:, [1 2 3 4 5]) - [1, 1, 1, 1, 0]) .* [ULM.SizeOfPixelZm, ULM.SizeOfPixelXm, ULM.SizeOfPixelZm, ULM.SizeOfPixelXm, 1]), TrajectoriesAll, 'UniformOutput', 0); % Store all trajectories into one single matrix
TrajectoriesHybrid = cellfun(@(x)((x(:, [1 2 3 4 5]) - [1, 1, 1, 1, 0]) .* [ULM.res, ULM.res, ULM.SizeOfPixelZm, ULM.SizeOfPixelXm, 1]), TrajectoriesAll, 'UniformOutput', 0); % Store all trajectories into one single matrix

% Build the different renderings based on PALA toolbox
MatrixTrack = []; MatrixVelz = []; MatrixVelnorm = [];

%% Tracked density matrix
MatrixTrack = ULM_Track2MatOut(TrajectoriesHybrid, ULM.SRsize(1:2) + [1 1] * 1, 'mode', '2D_tracks'); % build density matrix

if saveModeImg
    fprintf('\n Saving tracked density matrices MatrixTrack now please wait...\n');
    filenameImg = sprintf('MatrixTrack_%s.tif', timestamp);
    saveDataFigureULM(MatrixTrack, ULM, ...
        fullfile(saveFolder, 'Data', 'MatrixTrack.mat'), ...
        fullfile(saveFolder, 'Images', filenameImg), 'MatrixTrack');
    fprintf('\n Done \n');
end

%% Oriented velocity matrix
[MatrixTrack, MatrixVelz] = ULM_Track2MatOut(TrajectoriesHybrid, ULM.SRsize(1:2) + [1 1] * 1, 'mode', '2D_vel_z'); % build  matrix

if saveModeImg
    fprintf('\n Saving oriented velocity MatrixVelz now please wait...\n');
    filenameImg = sprintf('MatrixVelz_%s.tif', timestamp);
    Data.MatrixTrack = MatrixTrack; Data.MatrixVelz = MatrixVelz;
    saveDataFigureULM(Data, ULM, ...
        fullfile(saveFolder, 'Data', 'MatrixVelz.mat'), ...
        fullfile(saveFolder, 'Images', filenameImg), 'MatrixVelz');
    fprintf('\n Done \n');
end

%% Velocity norm matrix
[~, MatrixVelnorm] = ULM_Track2MatOut(TrajectoriesHybrid, ULM.SRsize + [1 1] * 1, 'mode', '2D_velnorm'); % build velocity norm matrix

if saveModeImg
    Data.MatrixTrack = MatrixTrack; Data.MatrixVelnorm = MatrixVelnorm;
    fprintf('\n Saving norm velocities filtered now please wait...\n'); drawnow;
    filenameImg = sprintf('MatrixVelnorm_filtered_%s.tif', timestamp);
    saveDataFigureULM(Data, ULM, ...
        fullfile(saveFolder, 'Data', 'MatrixVelnorm_filtered.mat'), ...
        fullfile(saveFolder, 'Images', filenameImg), 'MatrixVelnormfiltered');
    fprintf('\n Done \n');
end
