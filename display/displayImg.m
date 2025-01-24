function [ImgNorm, figNum] = displayImg(Img, varargin)
%% A function to display 2D images in a figure
    % Created by Baptiste Heiles on 2022/10/10
    % Last modified by Baptiste Heiles on 2025/01/03
    %
    % This function displays a 2D image in a new figure window. The image
    % can be normalized and displayed as a movie with various options.
    %
    % Usage:
    %   figNum = displayImg(Img) 
    %       Displays the image in a new figure with default settings.
    %       The image is log-normalized and displayed in a loop.
    %       The figure number is stored in figNum.BF.
    %
    %   figNum = displayImg(Img, normModes, dispModes)
    %       Displays the image with specified normalization and display modes.
    %       normModes: 'log', 'power', 'none'
    %       dispModes: 'loop', 'once'
    %       caxisval: 1x2 array specifying the color axis limits of the displayed image
    %
    % Inputs:
    %   Img: 2D array representing the image to be displayed
    %   normModes: (optional) string specifying the normalization mode
    %   dispModes: (optional) string specifying the display mode
    %   caxisval: (optional) 1x2 array specifying the color axis limits
    %
    % Outputs:
    %   ImgNorm: Normalized image
    %   figNum: Structure containing the figure number

    % Fetch and define variables based on input arguments
    switch nargin
        case 1
            normModes = 'log'; % Default normalization mode
            dispModes = 'loop'; % Default display mode
            caxisval = []; % Default color axis limits
        case 2
            normModes = varargin{1}; % User-specified normalization mode
            dispModes = 'loop'; % Default display mode
            caxisval = []; % Default color axis limits
        case 3
            normModes = varargin{1}; % User-specified normalization mode
            dispModes = varargin{2}; % User-specified display mode
            caxisval = []; % Default color axis limits
        case 4
            normModes = varargin{1}; % User-specified normalization mode
            dispModes = varargin{2}; % User-specified display mode
            caxisval = varargin{3}; % User-specified color axis limits
    end

    % Define or fetch a figNum variable to store the figure number
    f = gcf; % Get current figure
    figNum.BF = f.Number ; % Increment figure number

    % Create a new figure for displaying the image
    f = figure(figNum.BF); 
    f.InvertHardcopy = 'off'; % Preserve background color when printing
    f.Color = ones(1, 3); % Set figure background color to white
    tb = uicontrol(f, 'Style', 'togglebutton', 'String', 'Cancel'); % Add a cancel button
    drawnow; % Update figure window
    set(gca, 'Ydir', 'reverse'); % Reverse the y-axis direction
    assignin('base', 'figNum', figNum); % Assign figNum to the base workspace

    % Define the axes and time
    xaxis = 1:50:size(Img, 2); % Define x-axis ticks
    yaxis = 1:50:size(Img, 1); % Define y-axis ticks
    t_axis = []; % Initialize time axis

    % Pre-process Img
    Img = squeeze(Img); % Remove singleton dimensions

    % Normalize the image based on the specified normalization mode
    switch normModes
        case 'log'
            ImgNorm = 20 * log10(abs(Img) ./ max(abs(Img(:)))); % Log normalization
    
        case 'power'
            ImgNorm = (abs(Img) .^ 2) ./ mean(abs(Img) .^ 2, 2); % Power normalization
            ImgNorm = ImgNorm ./ max(ImgNorm(:));
        case 'none'
            ImgNorm = abs(Img); % No normalization
    end

    % Display the image
    switch dispModes
        case 'loop'
            loopTimeout = 0;
        case 'none'
            loopTimeout=19;
    end
    colorbar; colormap bone; % Set colormap to bone
    
    % Set axis labels and title
    xlabel('X-axis (pixels)');
    ylabel('Y-axis (pixels)');
    title('2D Image Display');

    while loopTimeout < 20  % Loop through the frames if needed
        axis image;
        set(gca, 'Ydir', 'reverse');
    
        for i_f = 1:size(ImgNorm, 3)
            if (get(tb, 'Value') == 1); display('Cancelled'); return; end
            imagesc(xaxis, yaxis, ImgNorm(:, :, i_f)); axis image
            switch normModes
                case 'log'
                    if isempty(caxisval)
                        caxis([-80 0]);
                    else
                        caxis(caxisval);
                    end
    
                case 'power'
                    if isempty(caxisval)
                        caxis([0 1]);
                    else
                        caxis(caxisval);
                    end
            end
    
            if ~isempty(t_axis)
                title(['t = ', num2str(i_f .* t_axis, '%3.3f'), ' s']);
            else
                title(['Frame = ', num2str(i_f)]);
            end
            colorbar
            drawnow; pause(0.001);

        end
    
        loopTimeout = loopTimeout + 1;
    end

end
