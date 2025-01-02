function o_dLos = Pixel2LoS_NoDistorsion(i_duvPixCoord, i_dKcam) %codegen
%% PROTOTYPE
% o_dLos = Pixel2LoS_NoDistorsion(i_duvPixCoord, i_dKcam) %codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function converting pixel coordinates in the image plane to line of sight
% unit vector for the given Pin-hole camera model (no distorsions).
% REFERENCE: Opromolla slides (file 2)
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% uv_pix: [2xN] arrays of N (u,v) pixels coordinates to convert
% Kcam: [3x3] Camera (non-dimensional) calibration matrix (intrinsic
%             parameters). Note: Focal length shall be normalized as
%             f/pixel_dim; the Camera centre shall be in pix.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% LoS: [3xN] Line of Sight vectors corresponding to (u, v) pixel
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-06-2023    Pietro Califano     Coded and verified.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Inclusion of model for distorsions

%% Function code

% Get number of pixels to convert
[nRows, Npixels] = size(i_duvPixCoord);
% uv_pix must be [2xNpixels]
uv_pix_array = i_duvPixCoord;

if isvector(i_dKcam)
    KcamMat = reshape(i_dKcam, 3, 3);
else
    KcamMat = i_dKcam;
end

if Npixels == 2
    uv_pix_array = i_duvPixCoord';
    [~, Npixels] = size(i_duvPixCoord);
elseif Npixels ~= 2 && nRows ~= 2
    error('uv_pix input not consistent')
end

% if distorted_flag == 0

% Retrieve Intrinsic camera parameters (no distorsions)
% Image centre position in (u,v) detector plane
cx = KcamMat(1, 3); % [pix]
cy = KcamMat(2, 3); % [pix]
% Normalized focal length (f/dim_pix)
f_nondim = KcamMat(1, 1); % [-]

% Static allocation
o_dLos = zeros(Npixels, 3);

% Evaluate LoS for every pixel
for i = 1:Npixels
    
    % Compute Normalized pixel components in image plane
    xNormPix = (uv_pix_array(1, i) - cx)/f_nondim;
    yNormPix = (uv_pix_array(2, i) - cy)/f_nondim;
    
    % Compute Bearing angles
    % Azimuth angle (positive rightward)
    Az = atan(xNormPix);
    cosAz = cos(Az);
    sinAz = xNormPix*cosAz;
    
    % Elevation angle (positive upward)
    El = atan(-yNormPix * cosAz);
    cosEl = cos(El);
    
    % Assuming Z as boresight
    o_dLos(i, :) = [cosEl*sinAz;
                   -sinAz
                   cosEl*cosAz];
end



end
