function Ap = postprocess_plateaus(A, E)
%postprocess_plateaus Replace upslope areas for plateaus with mean value 
%
%   Ap = postprocess_plateaus(A, E) "flattens" the upslope area computation for
%   plateaus in a DEM.  A is the upslope area for the DEM E, as computed by
%   upslope_area.  Ap is the same as A, except that the upslope area values for
%   each plateau are replace by the mean value for that plateau.
%
%   Example
%   -------
%
%       s = load('milford_ma_dem');
%       E = s.Zc;
%       R = dem_flow(E);
%       T = flow_matrix(E, R);
%       A = upslope_area(E, T);
%       Ap = postprocess_plateaus(A, E);
%       subplot(1,2,1)
%       imshow(log(A), [])
%       title('Upslope Area')
%       subplot(1,2,2)
%       imshow(log(Ap), [])
%       title('Upslope Area - Postprocessed')
%
%   See also upslope_area.

%   Steven L. Eddins
%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2007/08/06 19:46:57 $

no_downhill_neighbors = imerode(E, ones(3,3)) == E;

plateau_labels = bwlabel(no_downhill_neighbors);

s = regionprops(plateau_labels, 'PixelIdxList');

Ap = A;
for k = 1:numel(s)
    plateau_k_mean = mean(Ap(s(k).PixelIdxList));
    Ap(s(k).PixelIdxList) = plateau_k_mean;
end
