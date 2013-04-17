% computes p-field for anchor point.
% Inputs:
%     g: grid information (structure)
%     topo: initial topography (currently not used, but it may be used).
%     trend_map: linearly decreasing trend map (check trend_map_maker).
%     center_line_proxy_map: probability map decreasing from area of previous
% lobe simulated.
% Outputs:
%     pfield_start_point: Normalized p-field contructed considering all the 
% inputs. The map is normalized. 
% Written by Alejandro D. Leiva

function [pfield_startpoint]=make_pfield_for_startpoint_tau2(g,topo,trend_map,center_line_proxy_map)
lobe_prob_map_data=ones(g.ny,g.nx);

if size(center_line_proxy_map,1)==0;
    pfield_startpoint = trend_map./sum(sum(trend_map));
else

    %===========Base on topography===============================
    ele_max = max(max(topo));
    ele_min = min(min(topo));
    prob_map_ele = 1-(topo-ele_min)./(ele_max-ele_min);
      
    % normalize all the pmaps 
    prob_map_ele_norm = prob_map_ele./sum(sum(prob_map_ele));
    trend_map_norm = trend_map./sum(sum(trend_map));
    center_line_proxy_map_norm = center_line_proxy_map./sum(sum(center_line_proxy_map));
    lobe_prob_map_data_norm = lobe_prob_map_data./sum(sum(lobe_prob_map_data));

    tau_center_line = 4.0;
    tau_trend = 1.0;

    pA = 1/(g.nx*g.ny);   
    pADele = prob_map_ele_norm;   % These are pmaps - probablility of A occuring in each place assuming D, which is the elevation, migration, etc.
    pADtrend = trend_map_norm;    
    pADcenter = center_line_proxy_map_norm;
    pADdat = lobe_prob_map_data_norm;

    x0 = (1-pA)./pA;
    xcenter = zeros(g.ny, g.nx);
    xtrend = zeros(g.ny, g.nx);
    xall = zeros(g.ny, g.nx);
    prob_map_combined = zeros(g.ny, g.nx);

    for i=1:g.nx
        for j=1:g.ny
            if pADcenter(j,i)~=0&&pADtrend(j,i)~=0&&pADele(j,i)~=0&&pADdat(j,i)~=0
                xcenter(j,i)=(1-pADcenter(j,i))/pADcenter(j,i);
                xtrend(j,i)=(1-pADtrend(j,i))/pADtrend(j,i);
                %---------
                xall(j,i)=x0.*((xcenter(j,i)/x0)^tau_center_line)*((xtrend(j,i)/x0)^tau_trend);
                prob_map_combined(j,i)=1/(1+xall(j,i));
                %---------
            else    % If any of the probabilities =0, the whole thing =0
                prob_map_combined(j,i)=0;
            end
        end
    end

    sum_pmap_comb=sum(sum(prob_map_combined));
    pfield_startpoint=prob_map_combined/sum_pmap_comb;
end