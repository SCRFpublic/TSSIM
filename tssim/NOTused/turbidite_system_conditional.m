%Main script for surface-based turbidite lobe system simulation
%==========================================================================
%====== Main Program - Turbidite System Simulation (TSSIM)=================
%==========================================================================
%==========================================================================
%==== This program aims to simulate the stacking patterns observed in  ====
%==== distributary channel-lobe system using process-based statistics. ====
%==== It uses surface-based simulation, grid deformation, and Object-  ====
%==== based simulation.                                                ====

%==== by Alejandro Leiva, June 2009                                    ====
%==========================================================================
%==========================================================================
%==========================================================================

%--------------------------------------------------------------------------
%----------- 1. Clears and closes everything ------------------------------
%--------------------------------------------------------------------------
%clear;tic;clc;close all;
%--------------------------------------------------------------------------
%----------- 2. Looking for the path --------------------------------------
%--------------------------------------------------------------------------
dirname = get_dirname;
%--------------------------------------------------------------------------
%----------- 3. Simulation parameters -------------------------------------
%--------------------------------------------------------------------------
ifgu_threshold = 12000;
int_dep_rate = 5e-5;%Avg. dep. rate [m/yr] for time between lobe deposition 
nlobe = 6;
trend_map_range = 0.3;
s_source_windows = 60;
max_lobe_thickness_allowed = 3;%[m]
%----------- Flags --------------------------------------------------------
show_plots = false;
keep_plots_after_iteration = true;
consider_erosion = false; max_erosion_fraction = 0.2;
add_noise = false;
add_noise_ifgu = false;
smooth_surface = true; lmv.smoothing = 11;
conditional = true; n_wells = 3;
%--------------------------------------------------------------------------
%----------- 4. Loading initial surface------------------------------------
%--------------------------------------------------------------------------
%[base_topo]=flipud(read_earthvision_grid(path,g));
load base_topo;
[M N] = size(base_topo);
%--------------------------------------------------------------------------
%----------- 5. Making G-grid ---------------------------------------------
%--------------------------------------------------------------------------
verbose=0;
g.x0=0;
g.y0=0;
g.z0=0;
g.dx=2;
g.dy=2;
g.dz=1;
g.lz=1;
g.lx=N*g.dx;
g.ly=M*g.dy;
g = make_grid(g,0);
%path='Base_Surface_Model_2_600x500.txt';
%--------------------------------------------------------------------------
%----------- 6. Sediment source center ------------------------------------
%--------------------------------------------------------------------------
center_point = [round(g.nx/2) 1];
%--------------------------------------------------------------------------
%----------- 7. Get all of the CDF data from the process model ------------
%--------------------------------------------------------------------------
Input_path = strcat(dirname,'\code\');
lobe_length_cdf_data = [Input_path '\CDFS\lobe_length_cdf.txt'];
lobe_width_cdf_data = [Input_path '\CDFS\lobe_width_cdf.txt'];
max_lobe_thick_cdf_data = [Input_path '\CDFS\max_lobe_thick_cdf.txt'];
chan_detrend_totalL_to_width_cdf_data = [Input_path...
                              '\CDFS\channel_detrend_totL_over_W_cdf.txt'];
max_erosion_cdf_data = [Input_path '\CDFS\max_erosion_cdf.txt'];
intermediate_deposition_time_data = [Input_path...
                              '\CDFS\Intermediate_dep_time_cdf.txt'];
%----------- Making CDFs From of the files --------------------------------
lobe_length_cdf = read_cdf_data(lobe_length_cdf_data,verbose);
width_cdf = read_cdf_data(lobe_width_cdf_data,verbose);
max_lobe_thick_cdf = read_cdf_data(max_lobe_thick_cdf_data,verbose);
chan_detrend_totalL_to_width_cdf = read_cdf_data(...
                            chan_detrend_totalL_to_width_cdf_data,verbose);
max_erosion_cdf = read_cdf_data(max_erosion_cdf_data,verbose);
int_dep_time_cdf = read_cdf_data(intermediate_deposition_time_data,...
                    verbose);
%--------------------------------------------------------------------------
%---------- 8. Constructing trend map -------------------------------------
%--------------------------------------------------------------------------
trend_map_starting_point = trend_map_maker(g,center_point,trend_map_range);
[pfield_startpoint] = make_pfield_for_startpoint_tau2(g,base_topo,...
                      trend_map_starting_point,[]);
%--------------------------------------------------------------------------
%---------- 9. Memory Allocation ------------------------------------------
%--------------------------------------------------------------------------
sim_cat_cube_strata=zeros(g.ny,g.nx,nlobe*2);
sim_eros_cube_strata=zeros(g.ny,g.nx,nlobe*2);
sim_thick_cube_strata=zeros(g.ny,g.nx,nlobe*2);
surface_cube=zeros(g.ny,g.nx,nlobe*2);
tracklobe=-2*ones(1,nlobe*2);
%--------------------------------------------------------------------------
%---------- 10. Variables Initialization ----------------------------------
%--------------------------------------------------------------------------
tracklobe(1)= 1;%initial topography
layer_num = 2;
conditional_thickness = false;
core_type = 0;
lobe_num = 1; 
prev_sim_elev=base_topo;
surface_cube(:,:,1)=prev_sim_elev;
%--------------------------------------------------------------------------
%---------- 11. Reading hard data -----------------------------------------
%--------------------------------------------------------------------------
[conditional_interval hard_data_location warping_flag]= read_harddata(...
                                          conditional,n_wells,dirname,g);
%--------------------------------------------------------------------------
%---------- 12. Big loop --------------------------------------------------
%--------------------------------------------------------------------------
while lobe_num <= nlobe
    E = prev_sim_elev; 
%----------- Smoothening the surface --------------------------------------
    if smooth_surface == true
        E = lmv_surface(E,lmv.smoothing,g);
    end
%----------- Conditioning simulation for last lobe ------------------------
    if  conditional == true && nlobe - lobe_num + 1 <= ...
                           conditional_interval  && conditional_interval~=0
%----------- Gets data from file ------------------------------------------
        [length_core core_type type_contact well_last_data datai dataj ...
               pos_data] = get_core_data(hard_data_location,dirname,g);
%----------- Gets conditional simulation area -----------------------------           
        [map start_point checked_points I pos_i pos_j] = ...
                         get_simulation_area_conditional(E,[dataj...
                         datai],pfield_startpoint,g); 
        conditional_interval = conditional_interval -1;
        hard_data_location(2*M*N + pos_data) = hard_data_location...
                                                    (2*M*N + pos_data) + 1;
        conditional_thickness = true;
        if hard_data_location(2*M*N + pos_data) == well_last_data
            hard_data_location(pos_data) = 0;
        end
    else
%----------- Gets non conditional simulation area -------------------------
        [map start_point I pos_i pos_j] = get_simulation_area(E,...
                                        pfield_startpoint,show_plots,g,25);
        if  conditional == true
            wells_sim_area = sum(sum(map.*hard_data_location(:,:,1)~=0));
%----------- Checks if there is data in the SA ----------------------------
            if wells_sim_area ~= 0
                [length_core core_type type_contact well_last_data datai...
                         dataj pos_data] = get_core_data_in_SA(map, ...
                        hard_data_location,dirname,g);
                conditional_interval = conditional_interval -1;
                hard_data_location(2*M*N +pos_data) = hard_data_location...
                                                (2*M*N + pos_data) + 1;
                conditional_thickness = true;
                if hard_data_location(2*M*N + pos_data) == well_last_data
                    hard_data_location(pos_data) = 0;
                end
            end
        end
    end
    E = lmv_surface(map,11,g);           
    sim_area = map(:,:,1);
    contour_map = cont(sim_area);
%--------------------------------------------------------------------------
%---------- 13. Generating one lobe ---------------------------------------
%--------------------------------------------------------------------------

%----------- Gets allowable angle range -----------------------------------
    angles = get_angle(fliplr(start_point),contour_map,85);
    lobe_max_thick = interp1(max_lobe_thick_cdf(:,2), ...
                           max_lobe_thick_cdf(:,1),rand);
%----------- Checks is thickness is lower than max. thickness -------------                       
    if lobe_max_thick > max_lobe_thickness_allowed
        lobe_max_thick = max_lobe_thickness_allowed;
    end
    lobe_inside = 0;
%----------- Iterates until getting a lobe in the SA ----------------------                      
    while (lobe_inside==0)
        s_source_location = get_source_location(prev_sim_elev,[pos_j ...
                            pos_i],s_source_windows);
        [lobe_ind_map] = lobe_generator(start_point,...
            s_source_location,angles,g,verbose,Input_path);
        ind_lobe = length(find(lobe_ind_map)~=0);
        ind_lobe_inside = length(find(lobe_ind_map.*sim_area)~=0);
        ratio = ind_lobe_inside/ind_lobe;
        if ratio>0.8
            lobe_inside=1; 
            disp('Generating Object-based Lobe Simulation.... ')
            disp(strcat('Done with lobe number..',num2str(lobe_num),...
                    '!, there are..' ,num2str(nlobe-lobe_num),' left..'));
        end
    end
    sim_area_bi = sim_area>0;
    w = sim_area_bi+2*lobe_ind_map; 
    ww_sm = lmv_surface(w,31,g);
    if show_plots == true
        figure;imagesc(ww_sm);
    end
    lobe_ind_map = ww_sm > ww_sm(start_point(2),start_point(1));
%---------- Obtaining lobe channel ---------------------------------------- 
    s_source_location1(1) = g.nx-s_source_location(1);
    s_source_location1(2) = g.ny-s_source_location(2);
    start_point1(1) = g.nx-start_point(1);
    start_point1(2) = g.ny-start_point(2);
    d_chan = centerline_generator(s_source_location1,start_point1,g); 
    [across_dist] = find_nearest_point(d_chan,g,verbose);
    [chan_thick,chan_ind] = channel_shape_generator(across_dist,...
                 chan_detrend_totalL_to_width_cdf,lobe_max_thick,d_chan,g);
    ind_map =lobe_ind_map + fliplr(flipud(chan_ind));
%---------- Performs morphological closing --------------------------------
    se = strel('disk',7);
    closeBW = imclose(ind_map>0,se);
    sim_layer = closeBW;
%--------------------------------------------------------------------------    
%---------- 14. Lobe Thickeness Assignment --------------------------------
%--------------------------------------------------------------------------

%---------- if not conditional, noise can be added ------------------------
if conditional_thickness == false;
    thickness_map = lobe_thickness_assignment(closeBW,start_point,g);
%---------- Adds noise ----------------------------------------------------
    if add_noise == true
        noisy_data('hard_data_deposition.txt',thickness_map,0.5,g) 
        write_dimension_parfile_depo('sgsim_deposition.txt',g.nx,g.ny);
        dos('sgsim sgsim_deposition.txt');
        fid = fopen('sgs_deposition.txt','r');
        for i = 1:3 % skip the first 3 lines
            fgetl(fid);
        end
        cond_thickness_map = fscanf(fid,'%f');%concatenate data in 1D vect.
        cond_thickness_map = lobe_max_thick*reshape(cond_thickness_map,...
                             [g.nx  g.ny])';
        fclose(fid);
    else
        cond_thickness_map = lobe_max_thick*thickness_map;
    end
else
%---------- Conditional thickness assignment ------------------------------
    cond_thickness_map = lobe_cond_thickness_assignment(closeBW,...
                         start_point,dataj,datai,pos_data,length_core,...
                         max_lobe_thickness_allowed,g);
    conditional_thickness = false;
end
if show_plots == true
    figure;imagesc(cond_thickness_map);colorbar
end
%--------------------------------------------------------------------------
%--------- 15. Erosion ----------------------------------------------------
%--------------------------------------------------------------------------
    chan_e = interp1(max_erosion_cdf(:,2),max_erosion_cdf(:,1),rand);
%---------- Checking if max erosion is below threshold --------------------
    if chan_e > lobe_max_thick*max_erosion_fraction
        chan_e = lobe_max_thick*max_erosion_fraction;
    end
%---------- Computing erosion components ----------------------------------
    [profc,planc] = curvature(flipud(prev_sim_elev));  
    [G,ASP] = gradient8(flipud(prev_sim_elev));
    r_dir = river_direction(I,start_point);
    matrix_dir(:,:,1) = abs(ASP - r_dir(1));
    matrix_dir(:,:,2) = abs(ASP - r_dir(2));
    min_dir = min(matrix_dir,[], 3);
%---------- Calling functions the computes erosion map --------------------
    erosion_thick = erosion(sim_layer,profc,G,cond_thickness_map,...
                    min_dir,center_point,g,4/5);
    if show_plots == true
        figure;imagesc(erosion_thick);colorbar;title('Erosion Map');
    end
%---------- Adds noise in case it is required -----------------------------
    if add_noise == true && conditional == false
        noisy_data('hard_data_erosion.txt',erosion_thick,0.5,g) 
        write_dimension_parfile_erosion('sgsim_erosion.txt',g.nx,g.ny);
        dos('sgsim sgsim_erosion.txt');
        fid = fopen('sgs_erosion.txt','r');
        for i = 1:3
            fgetl(fid);
        end
        cond_erosion_map = fscanf(fid,'%f');  
        cond_erosion_map = reshape(cond_erosion_map,[g.nx  g.ny])';
        fclose(fid);
        cond_erosion_map = - chan_e*cond_erosion_map;
    else 
        cond_erosion_map =  - chan_e*erosion_thick;
    end
    if show_plots == true
        figure;imagesc(cond_erosion_map);
        colorbar;title('Conditional Erosion Map');
    end
    if consider_erosion == false
        cond_erosion_map =0;
    end
    tracklobe(layer_num)=lobe_num;
%---------- Making sure stacking rules are met ----------------------------
    dummy = surface_cube(:,:,layer_num-1)+cond_erosion_map;
    if layer_num>3
         dummy2 = surface_cube(:,:,layer_num-2);
         excp = dummy < dummy2;
         dummy(excp) = dummy2(excp);     
    end
%---------- Saving layers -------------------------------------------------
    surface_cube(:,:,layer_num-1) = dummy;
    surface_cube(:,:,layer_num) = surface_cube(:,:,layer_num-1)+...
                                  cond_thickness_map;
%--------------------------------------------------------------------------
%---------- 16. Fine Grain Unit -------------------------------------------
%--------------------------------------------------------------------------
    fines_dep_time = interp1(int_dep_time_cdf(:,2),...
                        int_dep_time_cdf(:,1),rand);
%---------- Computes layer thickness --------------------------------------
    fines_thick = int_dep_rate* fines_dep_time;%Avrg thickness of the unit
%---------- Checking if there is IFGU simulation, then simulates ----------
    if fines_dep_time > ifgu_threshold% Threshold for intermediate unit
        layer_num=layer_num+1;
        if add_noise_ifgu == true
            distribution = normrnd(fines_thick,0.2*fines_thick,100,1);
            fines_thick = fines_thick*ones(M,N);
            fid =  fopen('hist_ifgu.txt','w');
            fprintf(fid,'histogram\r\n');
            fprintf(fid,'%.0f\r\n',1);
            fprintf(fid,'histogram\r\n');
            for i = 1:100
                fprintf(fid,'%6.2f\r\n',distribution(i));
            end
            fclose(fid)
            noisy_data('hard_data_ifgu.txt',fines_thick,0.5,g)%writing data
            write_dimension_parfile_fgu('sgsim_ifgu.txt',g.nx,g.ny)
            dos('sgsim sgsim_ifgu.txt');
            fid = fopen('sgs_ifgu.txt','r');
            for i = 1:3
                fgetl(fid);
            end
            cond_ifgu_map = fscanf(fid,'%f'); 
            cond_ifgu_map = reshape(cond_ifgu_map,[g.nx  g.ny])';
            fclose(fid);
            cond_ifgu_map(cond_ifgu_map<prctile(cond_ifgu_map(:),25)) = 0;
        else 
            cond_ifgu_map = get_fgu(prev_sim_elev,fines_thick,0.5);
        end
        if show_plots == true
            figure;imagesc(cond_ifgu_map);colorbar;
            title('Intermediate Fine-grained Unit');
        end
%-------- Stacking layer (surface_cube) and assigning flag (tracklobe) ----        
        tracklobe(layer_num) = 0;    
        surface_cube(:,:,layer_num) = surface_cube(:,:,layer_num-1)+...
                                      cond_ifgu_map;
    end % end if fine-grained intermediate unit simulated
%--------------------------------------------------------------------------
%---------- 17. Re-starting variables -------------------------------------
%--------------------------------------------------------------------------
    prev_sim_elev = surface_cube(:,:,layer_num);
    layer_num = layer_num+1;
    lobe_num = lobe_num +1;
%---------- updating the pfield_startpoint --------------------------------
    s_source_location1(1) = g.nx-s_source_location(1);
    s_source_location1(2) = g.ny-s_source_location(2);
    start_point1(1) = g.nx-start_point(1);
    start_point1(2)= g.ny-start_point(2);
    d_chan = centerline_generator(s_source_location1,start_point1,g);
    for i=1:size(d_chan.ii,2)
        I(g.ny-d_chan.jj(i),g.nx-d_chan.ii(i))=1;
    end
    center_line_proxy_map = center_line_proxy(sim_layer);
    [pfield_startpoint] = make_pfield_for_startpoint_tau2(g,base_topo,...
                trend_map_starting_point,center_line_proxy_map);
    if show_plots == true
        figure;imagesc(pfield_startpoint);colorbar;
        title('P-field for Next Iteration');
    end
    if keep_plots_after_iteration == false
        close all;
    end
    clc;
end
%toc;
%---------- End of loop ---------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%---------- by Alejandro D. Leiva -----------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------