% Assigns conditional thickness to binary 2D geobody (closeBW).
% Inputs:
%     closeBW: Binary 2D geobody map.
%     start_point: Anchor point.
%     dataj: j-coord. data.
%     datai: i-coord. data.
%     pos_data: data positions (in a g.nx*g.ny matrix).
%     length_core: length of the conditional interval.
%     max_lobe_thickness_allowed: maximum lobe thickness allowed by user in 
% the simulation.
%     g: grid dimensions.
% Outputs:
%     cond_thickness_map: conditional map. 2D altitude map.
% Written by Alejandro D. Leiva, June '09
function cond_thickness_map = lobe_cond_thickness_assignment(closeBW,start_point,dataj,datai,pos_data,length_core,max_lobe_thickness_allowed,g)
cont_lobe = cont(closeBW);
[i_min j_min] = find_closest(cont_lobe,[dataj datai]);
closeBW_thick = lobe_thickness_assignment(closeBW,start_point,g);
closeBW_thick(closeBW_thick>0.95) = 1;
[i2_min j2_min] = find_closest(closeBW_thick,[j_min i_min]);
c_mult = 0.8; %convex multiplier for warping
while (true)
    originalMarks = [round(j_min*c_mult+j2_min*(1-c_mult)) round(i_min*c_mult+i2_min*(1-c_mult))];
    [closeBW1] = warpImage(closeBW, originalMarks, [dataj datai]);
    cond_thickness_map_mod = lobe_thickness_assignment(closeBW1,start_point,g);
    cond_thickness_map_mod = length_core*cond_thickness_map_mod/cond_thickness_map_mod(pos_data);
    if max(max(cond_thickness_map_mod)) < max_lobe_thickness_allowed
        break;
    else
        c_mult = c_mult - 0.05;
    end   
end
cond_thickness_map = cond_thickness_map_mod;