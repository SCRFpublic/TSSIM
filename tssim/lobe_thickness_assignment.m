% Assigns thickness to binary 2D geobody (closeBW). Returns thickness map
% Inputs:
%     closeBW: Binary 2D geobody map.
%     start_point: Anchor point.
%     g: grid dimensions.
% Outputs:
%     cond_map: conditional map. 2D altitude map from 0 to 1.
% Written by Alejandro D. Leiva, June '09
function thickness_map = lobe_thickness_assignment(closeBW,start_point,g)
dist_anchor = zeros(size(closeBW));
dist_anchor(start_point(2),start_point(1)) = 1;
dist_anchor = bwdist(dist_anchor);
dist_anchor = dist_anchor.*(closeBW>0);
dist_anchor(closeBW>0) =1-dist_anchor(closeBW>0)/max(max(dist_anchor)); %tunning thickness parameter (check this)
prox_dist = bwdist(1-closeBW);
thickness_map = prox_dist.*dist_anchor;
low_zones = find(thickness_map>0.01 & thickness_map<0.5);
thickness_map(low_zones) = thickness_map(low_zones)+0.2;
ww=zeros(g.ny,g.nx,2);
ww(:,:,1)=thickness_map;
ww_smooth=smooth3(ww,'box',[15 15 1]);
thickness_map=ww_smooth(:,:,1);
thickness_map = thickness_map/max(max(thickness_map));