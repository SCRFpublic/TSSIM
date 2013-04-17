function A = upslope_area(E, T)
%upslope_area Upslope area measurements for a DEM
%
%   A = upslope_area(E, T) computes the upslope area for each pixel of the
%   DEM matrix, E.  T is the sparse system of linear equations computed by
%   flow_matrix; it represents the distribution of flow from pixel to pixel.
%
%   A contains the upslope area for each corresponding pixel of E.
%
%   Note: Connected groups of NaN pixels touching the border are treated as
%   having no contribution to flow.
%
%   Reference: Tarboton, "A new method for the determination of flow
%   directions and upslope areas in grid digital elevation models," Water
%   Resources Research, vol. 33, no. 2, pages 309-319, February 1997. 
%
%   Algorithm notes: The Tarboton paper is not very specific about the
%   handling of plateaus.  For details of how plateaus are handled in this
%   code, see the algorithm notes for the function flow_matrix.  In
%   particular, see the subfunction plateau_flow_weights in flow_matrix.m.
%   
%   Example
%   -------
%
%       s = load('milford_ma_dem');
%       E = s.Zc;
%       R = dem_flow(E);
%       T = flow_matrix(E, R);
%       A = upslope_area(E, T);
%       imshow(log(A), [])
%
%   See also dem_flow, dependence_map, fill_sinks, flow_matrix,
%   postprocess_plateaus.

%   Steven L. Eddins
%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 2007/10/02 15:50:40 $

% Right-side vector is normally all ones, reflecting an equal contribution
% to water flow originating in each pixel.
rhs = ones(numel(E), 1);

% Connected groups of NaN pixels that touch the border do not contribute
% to water volume.
mask = border_nans(E);
rhs(mask(:)) = 0;

A = T \ rhs;
A = reshape(A, size(E));
