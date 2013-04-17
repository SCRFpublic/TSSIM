function [R, S] = D8_flow(E, i, j, d1)
%Compute downslope flow direction.

if nargin < 4
    d1 = 1;
end
[M, N] = size(E);

% Preprocess NaNs connected to the border.  Make them higher than the
% highest non-NaN value.
highest_value = max(E(:));
bump = min(1, eps(highest_value));
E(border_nans(E)) = highest_value + bump;

% Compute linear indices at desired locations.
e0_idx = (j - 1)*M + i;

% Table 1, page 311
% Row and column offsets corresponding to e1 and e2 for each
% table entry:
e1_row_offsets = [1 1 0 -1 -1 -1 0 1];
e1_col_offsets = [0 -1 -1 -1 0 1 1 1];


% Linear e1 offsets.
e1_linear_offsets = e1_col_offsets*M + e1_row_offsets;


% Initialize R and S values based on the first facet.
E0 = E(e0_idx);
E1 = E(e0_idx + e1_linear_offsets(1));

[R, S] = D8_pixel(E0, E1, d1);
R=6/8*2*pi;


for k = 2:8
    % Compute Rk and Sk corresponding to the k-th facet. Where Sk is positive
    % and greater than S, replace S and recompute R based on Rk.
    E1 = E(e0_idx + e1_linear_offsets(k));
    d1=[sqrt(2) 1 sqrt(2) 1 sqrt(2) 1 sqrt(2)];
    [Rk, Sk] = D8_pixel(E0, E1, d1(k-1));
    
    if Sk>S
        S=Sk;
        R=((7-k)/8)*2*pi;
    end
  
end



function [r, s] = D8_pixel(e0, e1, d1)
if nargin < 3
    d1 = 1;
end

s = (e0 - e1) / d1;           
r = 0;             






