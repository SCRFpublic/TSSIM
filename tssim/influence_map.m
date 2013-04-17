function I = influence_map(E, T, arg3, arg4)
%influence_map Influence map for pixel flow in a DEM
%   I = influence_map(E, T, i, j) calculates an influence map, I, for the DEM
%   matrix, E.  T is the pixel flow matrix as computed by flow_matrix.  i and
%   j are vectors containing the row and column coordinates of the starting
%   pixels.  Each element of the matrix I contains the amount of pixel flow
%   received from the starting pixels.
%
%   I = influence_map(E, T, BW) uses the nonzero pixels in the mask image BW
%   as the starting pixels for the calculation.
%
%   Note: Connected groups of NaN pixels touching the border are treated as
%   having no contribution to flow.
%
%   Reference: Tarboton, "A new method for the determination of flow
%   directions and upslope areas in grid digital elevation models," Water
%   Resources Research, vol. 33, no. 2, pages 309-319, February 1997. 
%
%   Example
%   -------
%
%       E = peaks;
%       R = dem_flow(E);
%       T = flow_matrix(E, R);
%       I = influence_map(E, T, 38, 28);
%       vis_map(I, E, 38, 28)
%
%   See also dependence_map, flow_matrix, upslope_area, vis_map.

%   Steven L. Eddins
%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2008/02/14 15:50:31 $



if nargin < 4
    BW = arg3;
    [i, j] = find(BW);

else
    i = arg3;
    j = arg4;
end

[M, N] = size(E);

rhs = zeros(numel(E), 1);
idx = (j-1)*M + i;
rhs(idx) = 1;

I = T \ rhs;
I = reshape(I, M, N);

    