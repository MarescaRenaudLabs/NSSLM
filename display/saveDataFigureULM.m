function results = saveDataFigureULM(Data, ULM, SavingPathData, SavingPathImg, mode)
    %% A function to compute and save ULM data from PALA format in a figure
    % Created by Baptiste Heiles on 2022/10/10
    % Last modified by Baptiste Heiles on 2025/01/23
    %
    % Data: The matrix or data structure containing ULM (Ultra Localizations Microscopy) or PALA results (e.g., localization density).
    % ULM: A structure with parameters describing the localization or super-resolution imaging process.
    % SavingPathData: The full file path (string) to save the numerical data.
    % SavingPathImg: The full file path (string) to save the figure image.
    % mode: A string specifying how to process the data (e.g., 'MatrixLoc').
    %
    % results: A structure or variable summarizing the saved data and any relevant metadata (e.g., paths or processing details).

    % Check if the directories exist, if not create them
    if ~isdir(fileparts([SavingPathData]))
        mkdir(fileparts([SavingPathData]));
    end

    if ~isdir(fileparts([SavingPathImg]))
        mkdir(fileparts([SavingPathImg]));
    end

    % Process the data and save it based on the specified mode
    switch mode
        case 'MatrixLoc' % in that case Data is localization
            % First save as a figure
            fdummy = figure(1);
            MatrixLoc = Data;
            SqMatrixLoc = sqrt(MatrixLoc);
            imagesc(SqMatrixLoc); adummy = get(fdummy, 'CurrentAxes'); axis(adummy, 'image');
            title('Localizations density squared'); caxis(adummy, [0 max(SqMatrixLoc(:))]);
            colormap(adummy, 'gray');
            clbarf = colorbar(adummy);
            clbarf.Label.String = 'Counts (A.U)';
            ScaleOfPixelm = [ULM.SizeOfPixelZm ULM.SizeOfPixelXm] .* ULM.res;
            adummy.DataAspectRatio(1:2) = ScaleOfPixelm(1:2); fdummy.Color = 'w'; fdummy.InvertHardcopy = 'off';
            adummy.XTickLabel = round(adummy.XTick .* ULM.SizeOfPixelXm ./ ULM.res .* 1E4) ./ 10;
            adummy.YTickLabel = round(adummy.YTick .* ULM.SizeOfPixelZm ./ ULM.res .* 1E4) ./ 10;
            adummy.XLabel.String = 'X (mm)';
            adummy.YLabel.String = 'Z (mm)';
            drawnow;
            fmdummy = getframe(fdummy);
            print(fdummy, '-dtiffn', [SavingPathImg, '.tif']);
            fdummy.InvertHardcopy = 'off';
            savefig(fdummy, [SavingPathImg(1:end - 4), '.fig'], 'compact');
            close(fdummy);
            % Second save as an image with integrated colorbar
            % integrate colorbar
            clbsize = [100, 10]; if clbsize > size(SqMatrixLoc) / 2, clbsize = [1 1]; end
            SqMatrixLoc(1:clbsize(1), 1:clbsize(2)) = repmat(linspace(max(SqMatrixLoc(:)), 0, clbsize(1))', 1, clbsize(2)); % add velocity colorbar
            WriteTif(SqMatrixLoc, gray(256), [SavingPathImg(1:end - 4), '_asImg.tif'], 'caxis', [0 max(SqMatrixLoc(:))]);
            % Third save Data
            save(SavingPathData, 'SqMatrixLoc', 'MatrixLoc');
        case 'MatrixTrack' % in that case Data is trajectories
            % First save as a figure
            fdummy = figure(1);
            SqMatrixTrack = sqrt(Data);
            imagesc(SqMatrixTrack); adummy = get(fdummy, 'CurrentAxes'); axis(adummy, 'image');
            title('Squared density of tracks'); caxis(adummy, [0 max(SqMatrixTrack(:))]);
            colormap(adummy, 'hot');
            clbarf = colorbar(adummy);
            clbarf.Label.String = 'Counts (A.U)';
            ScaleOfPixelm = [ULM.SizeOfPixelZm ULM.SizeOfPixelXm] .* ULM.res;
            adummy.DataAspectRatio(1:2) = ScaleOfPixelm(1:2); fdummy.Color = 'w'; fdummy.InvertHardcopy = 'off';
            adummy.XTickLabel = round(adummy.XTick .* ULM.SizeOfPixelXm ./ ULM.res .* 1E4) ./ 10;
            adummy.YTickLabel = round(adummy.YTick .* ULM.SizeOfPixelZm ./ ULM.res .* 1E4) ./ 10;
            adummy.XLabel.String = 'X (mm)';
            adummy.YLabel.String = 'Z (mm)';
            drawnow;
            fmdummy = getframe(fdummy);
            print(fdummy, '-dtiffn', [SavingPathImg, '.tif']);
            fdummy.InvertHardcopy = 'off';
            savefig(fdummy, [SavingPathImg(1:end - 4), '.fig'], 'compact');
            close(fdummy);
            % Second save as an image with integrated colorbar
            % integrate colorbar
            clbsize = [100, 10]; if clbsize > size(SqMatrixTrack) / 2, clbsize = [1 1]; end
            SqMatrixTrack(1:clbsize(1), 1:clbsize(2)) = repmat(linspace(max(SqMatrixTrack(:)), 0, clbsize(1))', 1, clbsize(2)); % add velocity colorbar
            WriteTif(SqMatrixTrack, hot(256), [SavingPathImg(1:end - 4), '_asImg.tif'], 'caxis', [0 max(SqMatrixTrack(:))]);
            % Third save Data
            MatrixTrack = Data;
            save(SavingPathData, 'SqMatrixTrack', 'MatrixTrack');
        case 'MatrixVelz' % in that case Data is trajectories
            % in that case Data is a struct with the two useful matrices
            MatrixTrack = Data.MatrixTrack; MatrixVelz = Data.MatrixVelz .* 1e3; % hardcoded conversion to mm/s
            % Filter
            MatrixVelzFilt = sqrt(MatrixTrack) .* sign(imgaussfilt(MatrixVelz, .8));
            velColormap = cat(1, flip(flip(hot(128), 1), 2), hot(128));
            % First save as figure
            fdummy = figure(1);
            imagesc(MatrixVelzFilt); adummy = get(fdummy, 'CurrentAxes'); axis(adummy, 'image');
            title('Filtered oriented velocity (mm/s)'); caxis(adummy, max(abs(MatrixVelzFilt(:))) * [-1 1]);
            colormap(adummy, velColormap);
            clbarf = colorbar(adummy);
            clbarf.Label.String = 'Velocity sq (mm/s)';
            ScaleOfPixelm = [ULM.SizeOfPixelZm ULM.SizeOfPixelXm] .* ULM.res;
            adummy.DataAspectRatio(1:2) = ScaleOfPixelm(1:2); fdummy.Color = 'w'; fdummy.InvertHardcopy = 'off';
            adummy.XTickLabel = round(adummy.XTick .* ULM.SizeOfPixelXm ./ ULM.res .* 1E4) ./ 10;
            adummy.YTickLabel = round(adummy.YTick .* ULM.SizeOfPixelZm ./ ULM.res .* 1E4) ./ 10;
            adummy.XLabel.String = 'X (mm)';
            adummy.YLabel.String = 'Z (mm)';
            drawnow;
            fmdummy = getframe(fdummy);
            print(fdummy, '-dtiffn', [SavingPathImg, '.tif']);
            fdummy.InvertHardcopy = 'off';
            savefig(fdummy, [SavingPathImg(1:end - 4), '.fig'], 'compact');
            close(fdummy);
            % Second save as an image with integrated colorbar
            % Place colorbars
            clbsize = [100, 10]; if clbsize > size(MatrixVelzFilt) / 2, clbsize = [1 1]; end
            MatrixVelzFilt(1:clbsize(1), 1:clbsize(2)) = repmat(linspace(-max(abs(MatrixVelzFilt(:))), max(abs(MatrixVelzFilt(:))), clbsize(1))', 1, clbsize(2)); % add velocity colorbar
            WriteTif(MatrixVelzFilt, velColormap, [SavingPathImg(1:end - 4), '_asImg.tif'], 'caxis', max(abs(MatrixVelzFilt(:))) * [-1 1]);

            % Third save data
            %                 save(SavingPathData,'MatrixVelzFilt','MatrixVelz');
        case 'MatrixVelnormfiltered' % in that case Data is trajectories
            % in that case Data is a struct with the two useful matrices
            MatrixTrack = Data.MatrixTrack; MatrixVelnorm = Data.MatrixVelnorm .* 1e3; % hardcoded conversion to mm/s
            % Filter
            vmax_disp = ceil(quantile(MatrixVelnorm(abs(MatrixVelnorm) > 0), .98) / 10) * 10;
            % Normalize matrix by vmax_disp
            MatrixVelnormfilt = MatrixVelnorm / vmax_disp; % normalization
            MatrixVelnormfilt = MatrixVelnormfilt .^ (1/1.5); MatrixVelnormfilt(MatrixVelnormfilt > 1) = 1;
            MatrixVelnormfilt = imgaussfilt(MatrixVelnormfilt, .5);
            MatrixVelnormfilt = ind2rgb(round(MatrixVelnormfilt * 256), jet(256)); % convert ind into RGB
            MatShadow = MatrixTrack; MatShadow = MatShadow ./ max(MatShadow(:) * .3); MatShadow(MatShadow > 1) = 1;
            MatrixVelnormfilt = MatrixVelnormfilt .* (MatShadow .^ .25);

            % First save as figure
            fdummy = figure(1);
            imagesc(MatrixVelnormfilt); adummy = get(fdummy, 'CurrentAxes'); axis(adummy, 'image');
            title('Filtered norm of velocity'); caxis(adummy, [0 max(MatrixVelnormfilt(:))]);
            colormap(adummy, jet(256));
            clbarf = colorbar(adummy);
            clbarf.Label.String = 'Velocity (AU)';
            ScaleOfPixelm = [ULM.SizeOfPixelZm ULM.SizeOfPixelXm] .* ULM.res;
            adummy.DataAspectRatio(1:2) = ScaleOfPixelm(1:2); fdummy.Color = 'w'; fdummy.InvertHardcopy = 'off';
            adummy.XTickLabel = round(adummy.XTick .* ULM.SizeOfPixelXm ./ ULM.res .* 1E4) ./ 10;
            adummy.YTickLabel = round(adummy.YTick .* ULM.SizeOfPixelZm ./ ULM.res .* 1E4) ./ 10;
            adummy.XLabel.String = 'X (mm)';
            adummy.YLabel.String = 'Z (mm)';
            drawnow;
            fmdummy = getframe(fdummy);
            print(fdummy, '-dtiffn', [SavingPathImg, '.tif']);
            fdummy.InvertHardcopy = 'off';
            savefig(fdummy, [SavingPathImg(1:end - 4), '.fig'], 'compact');
            close(fdummy);
            % Second save as an image with integrated colorbar
            % Place colorbar
            clbsize = [100, 10]; if clbsize > size(MatrixVelnorm) / 2, clbsize = [1 1]; end
            MatShadow(1:clbsize(1), 1:clbsize(2)) = repmat(linspace(0, 1, clbsize(2)), clbsize(1), 1);
            imwrite(MatrixVelnormfilt, [SavingPathImg(1:end - 4), '_asImg.tif']);
            % Third save data
            MatrixVelnorm = Data.MatrixVelnorm;
            save(SavingPathData, 'MatrixVelnormfilt', 'MatrixVelnorm');
        case 'MatrixVelnorm' % in that case Data is trajectories
            % in that case Data is a struct with the two useful matrices
            MatrixTrack = Data.MatrixTrack; MatrixVelnorm = Data.MatrixVelnorm .* 1e3; % hardcoded conversion to mm/s
            % Filter
            MatrixVelnormFilt = sqrt(MatrixVelnorm) .* sign(imgaussfilt(MatrixVelnorm, .8));
            velColormap = hot(256);
            % First save as figure
            fdummy = figure(1);
            imagesc(MatrixVelnormFilt); adummy = get(fdummy, 'CurrentAxes'); axis(adummy, 'image');
            title('Filtered norm velocity (mm/s)'); caxis(adummy, [0 max(MatrixVelnormFilt(:))]);
            colormap(adummy, velColormap);
            clbarf = colorbar(adummy);
            clbarf.Label.String = 'Velocity sq (mm/s)';
            ScaleOfPixelm = [ULM.SizeOfPixelZm ULM.SizeOfPixelXm] .* ULM.res;
            adummy.DataAspectRatio(1:2) = ScaleOfPixelm(1:2); fdummy.Color = 'w'; fdummy.InvertHardcopy = 'off';
            adummy.XTickLabel = round(adummy.XTick .* ULM.SizeOfPixelXm ./ ULM.res .* 1E4) ./ 10;
            adummy.YTickLabel = round(adummy.YTick .* ULM.SizeOfPixelZm ./ ULM.res .* 1E4) ./ 10;
            adummy.XLabel.String = 'X (mm)';
            adummy.YLabel.String = 'Z (mm)';
            drawnow;
            fmdummy = getframe(fdummy);
            print(fdummy, '-dtiffn', [SavingPathImg, '.tif']);
            fdummy.InvertHardcopy = 'off';
            savefig(fdummy, [SavingPathImg(1:end - 4), '.fig'], 'compact');
            close(fdummy);
            % Second save as an image with integrated colorbar
            % Place colorbars
            clbsize = [100, 10]; if clbsize > size(MatrixVelnormFilt) / 2, clbsize = [1 1]; end
            MatrixVelnormFilt(1:clbsize(1), 1:clbsize(2)) = repmat(linspace(0, max(abs(MatrixVelnormFilt(:))), clbsize(1))', 1, clbsize(2)); % add velocity colorbar
            WriteTif(MatrixVelnormFilt, velColormap, [SavingPathImg(1:end - 4), '_asImg.tif'], 'caxis', [0 max(MatrixVelnormFilt(:))]);

            % Third save data
            save(SavingPathData, 'MatrixVelnormFilt', 'MatrixVelnorm');
    end

end
