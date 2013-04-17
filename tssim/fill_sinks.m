function Ef = fill_sinks(E)
%fill_sinks Fill interior sinks in DEM
%
%   Ef = fill_sinks(E) fills interior sinks in a DEM.
%
%   Example
%   -------
%
%       s = load('milford_ma_dem');
%       E = s.Zc;
%       Ef = fill_sinks(E);
%       R = dem_flow(Ef);
%       T = flow_matrix(Ef, R);
%       A = upslope_area(Ef, T);
%       imshow(log(A), [])
%
%   See also imfill, upslope_area.

%   Steven L. Eddins
%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2007/08/03 14:49:58 $

Ef = imfill(E, 'holes');
