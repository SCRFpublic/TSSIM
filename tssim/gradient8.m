% 8-connected neighborhood gradient and aspect of a digital elevation model
%
% [G,ASP] = gradient8(DEM,cellsize)
%
% gradient8 returns the numerical steepest downward gradient and aspect 
% of a digital elevation model using an 8-connected neighborhood. 
%
% Input arguments
%       DEM       m x n matrix with elevation values
%       cellsize  horizontal spacing between data points (default = 1)
% 
% Output arguments
%       G         matrix with tangent of the gradient
%       ASP       matrix with aspect containing integers from 1 to 8
%                 according to the direction of the slope counted clockwise
%                 from top.
%
%                    8 1 2
%                    7 c 3
%                    6 5 4
%
%                 ASP is zero for cells without downward neighbor. 
%                  
% Example:
% 
% [X,Y,DEM] = peaks(100);
% res = abs(X(1,1)-X(1,2));
% [G,ASP] = gradient8(DEM,res);
% subplot(1,2,1)
% surf(X,Y,DEM,G); shading interp; camlight
% subplot(1,2,2)
% surf(X,Y,DEM,ASP); shading flat; camlight
%
% 
% Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% 

function [G,ASP] = gradient8(DEM,res)


if nargin == 1;
    res = 1;
else
    if ~isscalar(res) || res <= 0
        error('cellsize must be a positive scalar')
    end
end

siz    = size(DEM);

% pad dem with nans
DEM    = [nan(1,siz(2)+2); [nan(siz(1),1) DEM nan(siz(1),1)]; nan(1,siz(2)+2)];

% function handles for shifting DEM
neighfun = cell(8,2);

neighfun{1,1} = @(x,d) (x(2:end-1,2:end-1)-x(1:end-2,2:end-1))/d;
neighfun{2,1} = @(x,d) (x(2:end-1,2:end-1)-x(1:end-2,3:end))/hypot(d,d);
neighfun{3,1} = @(x,d) (x(2:end-1,2:end-1)-x(2:end-1,3:end))/d;
neighfun{4,1} = @(x,d) (x(2:end-1,2:end-1)-x(3:end,3:end))/hypot(d,d);

neighfun{5,1} = @(x,d) (x(2:end-1,2:end-1)-x(3:end,2:end-1))/d;
neighfun{6,1} = @(x,d) (x(2:end-1,2:end-1)-x(3:end,1:end-2))/hypot(d,d);
neighfun{7,1} = @(x,d) (x(2:end-1,2:end-1)-x(2:end-1,1:end-2))/d;
neighfun{8,1} = @(x,d) (x(2:end-1,2:end-1)-x(1:end-2,1:end-2))/hypot(d,d);

% preallocate arrays
G   = zeros(siz);
ASP = G;

% loop through neighbors
for neigh = (1:8);
    G2       = neighfun{neigh,1}(DEM,res);
    I        = G2>G;
    G(I)     = G2(I);
    ASP(I)   = neigh;
end

