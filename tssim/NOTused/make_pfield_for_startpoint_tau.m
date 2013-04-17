function [pfield_startpoint]=make_pfield_for_startpoint_tau(base_topo,prob_map_mg, prob_map_pg,lobe_prob_map_data,sed_source_loc,g,topo_flag, condition_flag)

%This function combines migration, progradation, topographic (if desired),
%and data (if conditional) CDFs to get a probability map for anchor point
%locations.

% In this version, each new point is independent of the
% previously-simulated point.

% Written by: Holly Michael
% Date: 2008 Summer


%first calculate the probability map based on current elevation
% the lower the elevation, the higher the probability. This
% can be turned on or off by the topo_flag.
%===========Base on topography===============================
ele_max=max(max(base_topo));
ele_min=min(min(base_topo));
prob_map_ele=1-(base_topo-ele_min)./(ele_max-ele_min); % Added the 0.0001 at the end to avoid a divide by zero error when the initial surface is zero.
      
% normalize all the pmaps first (?)
prob_map_ele_norm=prob_map_ele./sum(sum(prob_map_ele));
prob_map_mg_norm=prob_map_mg./sum(sum(prob_map_mg));
prob_map_pg_norm=prob_map_pg./sum(sum(prob_map_pg));
lobe_prob_map_data_norm=lobe_prob_map_data./sum(sum(lobe_prob_map_data));

sum(sum(prob_map_ele_norm))
sum(sum(prob_map_mg_norm))
sum(sum(prob_map_pg_norm))
sum(sum(lobe_prob_map_data_norm))

tau_mg=1.0;
tau_pg=1.0;
tau_ele=2.0;
tau_dat=1.0;

pA=1/(g.nx*g.ny);   % Think this is supposed to be the probability of the startpoint occurring in a given cell. Check this.
pADele=prob_map_ele_norm;   % These are pmaps - probablility of A occuring in each place assuming D, which is the elevation, migration, etc.
pADmg=prob_map_mg_norm;    
pADpg=prob_map_pg_norm;
pADdat=lobe_prob_map_data_norm;

x0=(1-pA)./pA;
xmg=zeros(g.ny, g.nx);
xpg=zeros(g.ny, g.nx);
xele=zeros(g.ny, g.nx);
xdat=zeros(g.ny, g.nx);

for i=1:g.nx
    for j=1:g.ny
        if pADmg(j,i)~=0&pADpg(j,i)~=0&pADele(j,i)~=0&pADdat(j,i)~=0
            xmg(j,i)=(1-pADmg(j,i))/pADmg(j,i);
            xpg(j,i)=(1-pADpg(j,i))/pADpg(j,i);

            if topo_flag==1
                xele(j,i)=(1-pADele(j,i))/pADele(j,i);
            end

            if condition_flag==1
                xdat(j,i)=(1-pADdat(j,i))/pADdat(j,i);
            end
            %---------
            xall(j,i)=x0.*((xmg(j,i)/x0)^tau_mg)*((xpg(j,i)/x0)^tau_pg);
            if topo_flag==1
                xall(j,i)=xall(j,i)*((xele(j,i)/x0)^tau_ele);
            end

            if condition_flag==1
                xall(j,i)=xall(j,i)*((xdat(j,i)/x0)^tau_dat);
            end
            prob_map_combined(j,i)=1/(1+xall(j,i));
            %---------
        else    % If any of the probabilities =0, the whole thing =0
            prob_map_combined(j,i)=0;
        end
    end
end

sum_pmap_comb=sum(sum(prob_map_combined));
pfield_startpoint=prob_map_combined/sum_pmap_comb;


% figure;
% imagesc(g.x,g.y,xdat);
% title('Xdat');
% colorbar;
% xlabel('x [m]');
% ylabel('y [m]');
% axis xy;
% axis equal
% figure;
% imagesc(g.x,g.y,xele);
% title('Xele');
% colorbar;
% xlabel('x [m]');
% ylabel('y [m]');
% axis xy;
% axis equal
% figure;
% imagesc(g.x,g.y,xpg);
% title('Xpg');
% colorbar;
% xlabel('x [m]');
% ylabel('y [m]');
% axis xy;
% axis equal
% figure;
% imagesc(g.x,g.y,xmg);
% title('Xmg');
% colorbar;
% xlabel('x [m]');
% ylabel('y [m]');
% axis xy;
% axis equal

% figure;
% imagesc(g.x,g.y,prob_map_mg_norm);
% title('Normalized Migration Probability Map');
% colorbar;
% xlabel('x [m]');
% ylabel('y [m]');
% axis xy;
% axis equal
% 
% figure;
% imagesc(g.x,g.y,prob_map_pg_norm);
% title('Normalized Progradation Probability Map');
% colorbar;
% xlabel('x [m]');
% ylabel('y [m]');
% axis xy;
% axis equal
% 
% figure;
% imagesc(g.x,g.y,prob_map_ele_norm);
% title('Normalized Elevation Probability Map');
% colorbar;
% xlabel('x [m]');
% ylabel('y [m]');
% axis xy;
% axis equal
% 
% figure;
% imagesc(g.x,g.y,lobe_prob_map_data_norm);
% title('Normalized Data Probability Map');
% colorbar;
% xlabel('x [m]');
% ylabel('y [m]');
% axis xy;
% axis equal

% figure;
% imagesc(g.x,g.y,pfield_startpoint);
% title('Normalized Combined Probability Map');
% colorbar;
% xlabel('x [m]');
% ylabel('y [m]');
% axis xy;
% axis equal

