function vis_dem_flow(E, R, S)
%vis_dem_flow Visual pixel flow directions and magnitude on a DEM
%
%   vis_dem_flow(E, R, S) displays a DEM (E) as a grayscale image and
%   superimposes a quiver plot showing the direction (R) and magnitude (S) of
%   pixel flow directions.
%
%   Note: When R, S, and E contain more than 50 rows, vis_dem_flow crudely
%   downsamples R and S to avoid trying to display too many quiver arrows
%   simultaneously.
%
%   Example
%   -------
%
%       E = peaks;
%       [R, S] = dem_flow(E);
%       vis_dem_flow(E, R, S);
%
%   See also dem_flow.

%   Steven L. Eddins
%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2007/08/02 20:59:13 $


[M, N] = size(E);

imshow(E, [], 'InitialMagnification', 'fit')
hold on

[x, y] = meshgrid(1:N, 1:M);

if M > 50
    % Subsample S and R so we only try to plot no more than about 50 quiver
    % arrows vertically. 
    delta = M/50;
    y = y(round(1:delta:end), round(1:delta:end));
    x = x(round(1:delta:end), round(1:delta:end));
    S = S(round(1:delta:end), round(1:delta:end));
    R = R(round(1:delta:end), round(1:delta:end));
end

quiver(x, y, S.*cos(R), -S.*sin(R), 2, 'y')
hold off
