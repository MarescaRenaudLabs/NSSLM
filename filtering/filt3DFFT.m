function [filteredImg, fft_Img_filt, theta_fandef] = filt3DFFT(Img, filt_v, units, signkz, signkx)
    %% Function to filter in the time space domain
    % Created by Baptiste Heiles on 2024/04/24
    % Last modified by Baptiste Heiles on 2024/05/13
    %
    % This function will apply a 3D FFT and filter any moving objects according
    % to their velocity signature in the 3D Fourier space.
    %
    % Img: 3D array containing images with space x space x time
    % filt_v: 2x1 array containing filter boundaries in units same as input
    %    variable unit (m/s)
    % units: 3x1 array containing pixel units of image Img
    % signkz: int containing sign of filter along rows (z direction), -1
    %     towards top of the image, 1 towards the bottom, 0 for all direction
    % signkx: int containing sign of filter along columns (x direction), -1
    %     towards top of the image, 1 towards the bottom, 0 for all direction
    %
    % filteredImg: matrix containing filtered image
    % fft_Img_filt: array containing 3D FFT of input Img
    % theta_fandef: array with dimensions matching the input Img containing the
    %    mask used to filter velocities based on input

    fft_Img_filt = fftn(Img, size(Img));

    % Fetch speed limits
    komega_bound = filt_v; % m/s speed limits

    % Create Fourier space matching Img size
    nx = size(fft_Img_filt, 2);
    fx = 1 ./ units(2); % Spatial frequency in x direction
    nz = size(fft_Img_filt, 1);
    fz = 1 ./ units(1); % Spatial frequency in z direction
    nt = ceil(size(fft_Img_filt, 3));
    ft = 1 ./ units(3); % Temporal frequency
    [kx, kz, omega] = meshgrid(linspace(-fx ./ 2, fx / 2, nx), linspace(-fz ./ 2, fz / 2, nz), linspace(-ft ./ 2, ft / 2, nt));
    rho = sqrt(kz .^ 2 + kx .^ 2); % Define radius in Fourier space

    % Create mask to perform filtering
    theta_fandef = ones(size(fft_Img_filt)); % Initialize mask with ones

    theta_fandef(sign(omega .* kz) == signkz) = 0; % filter based on orientation along rows
    theta_fandef(sign(omega .* kx) == signkx) = 0; % filter based on orientation along columns
    theta_fandef(omega .* kz == 0) = 1; % Retain zero frequency components in kz
    theta_fandef(omega .* kx == 0) = 1; % Retain zero frequency components in kx

    % Now filter based on the input speeds
    theta_fandef(abs(omega) ./ rho > komega_bound(2)) = 0; % Inside border of the disk, larger speeds
    theta_fandef(abs(omega) ./ rho <= komega_bound(1)) = 0; % Outisde border of the disk, smaller speeds

    % Hamming filter on the fan to prevent artifacts
    win_z = round(0.1 .* nz); win_x = round(0.1 .* nx); win_t = round(0.01 .* nt);
    out = bsxfun(@times, bsxfun(@times, hamming(win_z), hamming(win_x).'), permute(hamming(win_t), [3 2 1]));
    theta_fandef = imfilter(theta_fandef, out, 'replicate');

    % Inverse transform taking account fft-shift of the fftn in Matlab (our filter
    %     is designed assuming kz=kx=omega=0 at center of the matrix)
    filteredImg = ifftn(fftshift(fft_Img_filt) .* theta_fandef);

end
