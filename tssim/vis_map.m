function vis_map(M, E, arg3, arg4)
%vis_map Visualize influence or dependence map on a DEM
%
%   vis_map(M, E, BW) displays an influence or dependence map (M) superimposed
%   on a DEM image (E).  BW is a binary image specifying the set of pixels from
%   which to measure influence or dependence.  M can be an influence map as
%   computed by influence_map, or it can be a dependence map as computed by
%   dependence_map.  Starting or ending pixel locations are shown in blue, and
%   pixels with nonzero influence or dependence values are shown in transparent
%   green.
%
%   vis_map(M, E, i, j) uses vectors i and j as row and column coordinates
%   for the starting pixel locations.
%
%   Example
%   -------
%
%       E = peaks;
%       R = dem_flow(E);
%       T = flow_matrix(E, R);
%       D = dependence_map(E, T, 12, 24);
%       vis_map(D, E, 12, 24)
%       title('Dependence map')
%
%   See also dependence_map, influence_map, upslope_area.

%   Steven L. Eddins
%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2007/08/02 20:59:13 $

if nargin < 4
    BW = arg3;
    
else
    i = arg3;
    j = arg4;
    BW = false(size(M));
    BW((j - 1)*size(M,1) + i) = true;
end

% Make and display an RGB image that is the autoscaled E matrix with the
% starting pixel locations shown in a shade of blue.
source_color = im2uint8([.5 .5 1]);
E = im2uint8(mat2gray(E));
red = E;
green = E;
blue = E;
red(BW) = source_color(1);
green(BW) = source_color(2);
blue(BW) = source_color(3);
rgb_bottom = cat(3, red, green, blue);
imshow(rgb_bottom, 'InitialMag', 'fit')

% Make a second RGB image that is a constant green.
rgb_top = zeros(size(M,1), size(M,2), 3, 'uint8');
rgb_top(:,:,2) = 255;

% Turn the influence or dependence map into an AlphaData channel to be used
% to display with the green image.
M(BW) = 0;
M = imadjust(mat2gray(M), [0 1], [0 .6], 0.5);

image('CData', rgb_top, 'AlphaData', M);
hold off
