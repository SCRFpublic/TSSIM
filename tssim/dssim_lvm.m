function [dssim_thick_map]=dssim_lvm(smooth_sim_thick,next_well_data,sgems_code_path,sgems_output_path_parfile,sgems_output_path_normal,matlab_output_path,matlab_code_path,g2s_code_path,layer_num,layer_flag,sillVal,cond, g)

% This program does direct sequential simulation (via SGEMS) using simple kriging
% with a locally-varying mean to get a smoothed, data matching map of
% thicknesses. 

%Input: 
% smooth_sim_thick: map of already smoothed (with a moving average) categorical thicknesses
% next_well_data: num_wells x 5 array of Well ID, X, Y, facies, max
    % thickness, where max thickness is the thickness of the whole facies unit
    % in the well log plus the max erosion thickness, OR the max facies
    % thickness in the TIs plus max erosion thickness (output from
    % get_next_data_layer
% sgems_code_path: path to the sgems multi ti code
% sgems_output_path_parfile: path to the sgems output folder, with slashes backwards from normal
% sgems_output_path_normal: path to the sgems output folder, with normal slashes
% matlab_output_path: matlab output path
% matlab_code_path: matlab code path
% g2s_code_path: gslib to binary conversion code path
% layer_num: layer number of simulation
% layer_flag: 0 if intermediate unit, 1 if lobe, 2 if erosion
% sillVal: value of the sill of the Gaussian variogram. Should be lower for erosion probably.
% cond: 1 if conditional simulation, 0 if not
% g: grid

% Outputs:
% dssim_thick_map: map of dssim-simulated thicknesses, based on thickness
% categories and using data where available. 

%---------------------------------------------------
% Create data file 

    hard_data_filename_gslib='hard_data_dssim.dat';
    hard_data_filename_binary='hard_data_dssim.out';

    fid=fopen(strcat(matlab_output_path,hard_data_filename_gslib),'w');
    if(fid==-1),
        fprintf('write_xydat_kriging_grid: Error opening file %s\n',fname);
        error('cannot open file');
    end
    fprintf(fid,'Description: Thickness Hard Data\n');
    fprintf(fid,'%.0f\n',3);
    fprintf(fid,'X\n');
    fprintf(fid,'Y\n');
    fprintf(fid,'Thickness\n');
    
    if cond==1  % If it's a conditional simulation
        % if it's a lobe layer:
        if layer_flag==1
            for wnum=1:size(next_well_data,1)
                if next_well_data(wnum,5)~=999  % if there's an interval left in the well
                    if next_well_data(wnum,7)~=0&next_well_data(wnum,8)~=-1&next_well_data(wnum,8)~=-3  % if it's not bypass
                        fprintf(fid,'%.0f %.0f %8.6f\n',next_well_data(wnum,2)-1,next_well_data(wnum,3)-1,next_well_data(wnum,5));
                    else    % if it is bypass
                        fprintf(fid,'%.0f %.0f %8.6f\n',next_well_data(wnum,2)-1,next_well_data(wnum,3)-1,0);
                    end
                else % if there's no data left in the well
                    fprintf(fid,'%.0f %.0f %8.6f\n',next_well_data(wnum,2)-1,next_well_data(wnum,3)-1,0);
                end
            end
            %             % also set the far edge to zero for lobe deposition:
            %             for cell_num=1:g.nx
            %                 fprintf(fid,'%.0f %.0f %8.6f\n',cell_num-1,1-1,0);
            %             end

            % if it's an intermediate finegrained unit:
        elseif layer_flag==0

            for wnum=1:size(next_well_data,1)
                if next_well_data(wnum,5)~=999  % if there's an interval left in the well
                    if next_well_data(wnum,4)==0&next_well_data(wnum,7)==0  % if it's finegrained and bypass
                        fprintf(fid,'%.0f %.0f %8.6f\n',next_well_data(wnum,2)-1,next_well_data(wnum,3)-1,0);
                    elseif next_well_data(wnum,4)==1  % if it's coarse-grained
                        fprintf(fid,'%.0f %.0f %8.6f\n',next_well_data(wnum,2)-1,next_well_data(wnum,3)-1,0);
                    else    % if it is fine-grained and not bypass
                        fprintf(fid,'%.0f %.0f %8.6f\n',next_well_data(wnum,2)-1,next_well_data(wnum,3)-1,next_well_data(wnum,5));
                    end
                else % if there's no data left in the well
                    fprintf(fid,'%.0f %.0f %8.6f\n',next_well_data(wnum,2)-1,next_well_data(wnum,3)-1,0);
                end
            end
           
        elseif layer_flag==2    % If it's erosion
            % all data points are designated as zero erosion except those that
            % have positive erosion thicknesses. erosion_depth_set, which is
            % read into next_well_data column 6, is 0 for coarse-grained with
            % no erosion, positive for coarse-grained with erosion (erosion
            % depth), and -1 for a fine-grained unit. Multiply by -1 to get the
            % actual erosion depth.
            for wnum=1:size(next_well_data,1)
                if next_well_data(wnum,6)~=999  % if there's an interval left in the well
                    if next_well_data(wnum,6)>0
                        fprintf(fid,'%.0f %.0f %8.6f\n',next_well_data(wnum,2)-1,next_well_data(wnum,3)-1,(-1*next_well_data(wnum,6)));
                    else    % if no erosion specified, set as 0
                        fprintf(fid,'%.0f %.0f %8.6f\n',next_well_data(wnum,2)-1,next_well_data(wnum,3)-1,0);
                    end
                else % if there's no data left in the well
                    fprintf(fid,'%.0f %.0f %8.6f\n',next_well_data(wnum,2)-1,next_well_data(wnum,3)-1,0);
                end
            end
            
             % because the data here is often all zero and we have to have a
            % nonzero data point (SGEMS requirement), assign some erosion at the sediment
            % source.
            %sed_source=[round(g.nx/2) g.ny];
            %fprintf(fid,'%.0f %.0f %8.6f\n',sed_source(1)-1,sed_source(2)-1,.1);

        end % end choosing which type of layer it is
    end % end if conditional
    % In all cases, set the far edge to something. This helps it not
    % have bumps at the edges.
    if layer_flag==1|layer_flag==0
        for cell_num=1:g.nx
            fprintf(fid,'%.0f %.0f %8.6f\n',cell_num-1,1-1,.02);
        end
    elseif layer_flag==2    % erosion
        for cell_num=1:g.nx
            fprintf(fid,'%.0f %.0f %8.6f\n',cell_num-1,1-1,0);
        end
    end

    fclose(fid);
    convert_to_binary('geoeas2sgems.par',strcat(matlab_output_path,hard_data_filename_gslib),strcat(matlab_output_path,hard_data_filename_binary),g2s_code_path,matlab_code_path,0,'hard_grid_dssim',g);

% Get parameters for DSSIM. The sill of the variogram should be equal to
% the variance of the histogram used, which will be the locally-varying
% mean  map. Actually, that's not the case, because the variance should be
% on the deviations from the LVM map, not the map itself. Since we don't
% have this information, let's take it as 1/10 the variance of the map.
% This can also be done as an input, but for now do it this way.

cat_vect=[];    % This is the LVM map put into vector format so that we can take the variance.
% first make the map into a vector:
for row_num=1:g.ny
    cat_vect=horzcat(cat_vect,smooth_sim_thick(row_num,:));
end
%size(cat_vect)
%g.nx*g.ny
sill=var(cat_vect)/10

% If decide to do it as a separate input:
%sill=sillVal;


% Also need the max and min of the lvm and the data for tail extrapolation
max_lvm=max(max(smooth_sim_thick));
min_lvm=min(min(smooth_sim_thick));
if cond==1
    if layer_flag==0 | layer_flag==1
        max_data=max(next_well_data(:,5));
        min_data=min(next_well_data(:,5));
    elseif layer_flag==2
        max_data=max(next_well_data(:,6));
        min_data=min(next_well_data(:,6));
    end

    % The value to use for the max in the simulation
    max_max=max(max_lvm,max_data);
    max_use=max_max+abs(max_max/20);
    min_min=min(min_lvm,min_data);
    min_use=min_min-abs(min_min/20);
else
    max_use=max_lvm+abs(max_lvm/20);
    min_use=min_lvm-abs(min_lvm/20);
end

% convert the smoothed thickness map to binary:
% first change X,Y to cell coords, not actual locations:
%----------
% get the locally-varying mean map from the smoothed simulated thickness
% category map. (Note that should change these thicknesses to the average
% value in each category, not the max!)
smooth_sim_thick_filename_gslib='smooth_sim_thick_file.dat';
smooth_sim_thick_filename_binary='smooth_sim_thick_file.out';
fid=fopen(strcat(matlab_output_path,smooth_sim_thick_filename_gslib),'w');
if(fid==-1)
    fprintf('write_xydat_grid: Error opening file %s\n',fname);
    error('cannot open file');
end

fprintf(fid,'Description: Thickness map: average thickness in each category\n');

fprintf(fid,'%d\n',3);
fprintf(fid,'Field_1:_x\n');
fprintf(fid,'Field_2:_y\n');
fprintf(fid,'Field_3:_thickness_or_lobe\n');

for j=1:g.ny
    for i=1:g.nx
        fprintf(fid,'%.6f %.6f %.6f\n',i-1,j-1,smooth_sim_thick(j,i));
    end
end
fclose(fid);

%-----------
convert_to_binary('geoeas2sgems.par',strcat(matlab_output_path,smooth_sim_thick_filename_gslib),strcat(matlab_output_path,smooth_sim_thick_filename_binary),g2s_code_path,matlab_code_path,0,'smooth_sim_thick',g);

binary_smooth_thicknesses_file=sprintf('smooth_sim_thick::');
hard_data_grid=sprintf('hard_grid_dssim"');
%---------------------------------------------------
% Now write the log file for running sgems   
parfile_name='dssim_parfile.log';
sim_name=strcat('dssim_data_map_',num2str(layer_num));
sim_filename=strcat(sim_name,'.dat');

% Took forever to get this to work. Tabs work differently each time. Seems
    % best to create variables above, using the \t followed by the next
    % character in the command line string.

% Write the parameter log file:
fid2=fopen(strcat(sgems_code_path,parfile_name),'wt');

if(fid2==-1),
    fprintf('write_xydat_logfile: Error opening file %s\n',parfile_name);
    error('cannot open file');
end

fprintf(fid2,'%s%s%d%s%d%s\n','NewCartesianGrid ','Sim_Grid::',g.nx,'::',g.ny,'::1::1::1::1.0::0::0::0');
fprintf(fid2,'%s %s\n','LoadObjectFromFile ',strcat(matlab_output_path,smooth_sim_thick_filename_binary,'::All'));
fprintf(fid2,'%s %s\n','CopyProperty ',strcat(binary_smooth_thicknesses_file,'Field_3:_thickness_or_lobe::Sim_Grid::Thickness::0::0'));

datafile_name=sprintf('hard_data_dssim');
fprintf(fid2,'%s %s\n','LoadObjectFromFile ',strcat(matlab_output_path,hard_data_filename_binary,'::All'));

if cond==1
%     datafile_name=sprintf('hard_data_dssim');
%     fprintf(fid2,'%s %s\n','LoadObjectFromFile ',strcat(matlab_output_path,hard_data_filename_binary,'::All'));

    if layer_flag==1 |layer_flag==0
        fprintf(fid2,'%s %s\n','RunGeostatAlgorithm ',['dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="Sim_Grid"  />     <Property_Name  value="',sim_name,'" />     <Nb_Realizations  value="1" />     <Seed  value="',num2str(ceil(100000+rand*100000)),'" />     <Kriging_Type  value="Kriging with Localy Varying Mean (LVM)"  />     <Local_Mean_Property  value="Thickness"  />     <Hard_Data  grid="',hard_data_grid,'   property="Thickness"  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="10" />     <Search_Ellipsoid  value="70 70 1  0 0 0" />    <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="0" filename =""   grid ="Sim_Grid"  property ="Thickness">  <LTI_type  function ="Power"  extreme ="0"  omega ="3" />  <UTI_type  function ="Power"  extreme ="',num2str(max_use),'"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="1e-006" structures_count="1"  >    <structure_1  contribution="',num2str(sill),'"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>']);
    elseif layer_flag==2
        fprintf(fid2,'%s %s\n','RunGeostatAlgorithm ',['dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="Sim_Grid"  />     <Property_Name  value="',sim_name,'" />     <Nb_Realizations  value="1" />     <Seed  value="',num2str(ceil(100000+rand*100000)),'" />     <Kriging_Type  value="Kriging with Localy Varying Mean (LVM)"  />     <Local_Mean_Property  value="Thickness"  />     <Hard_Data  grid="',hard_data_grid,'   property="Thickness"  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="10" />     <Search_Ellipsoid  value="70 70 1  0 0 0" />    <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="0" filename =""   grid ="Sim_Grid"  property ="Thickness">  <LTI_type  function ="Power"  extreme ="',num2str(min_use),'"  omega ="3" />  <UTI_type  function ="No extrapolation"  extreme ="0"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="1e-006" structures_count="1"  >    <structure_1  contribution="',num2str(sill),'"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>']);
   %                               RunGeostatAlgorithm  dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="SST_Erosion"  />     <Property_Name  value="test_lvm_1" />     <Nb_Realizations  value="1" />     <Seed  value="14071793" />     <Kriging_Type  value="Kriging with Localy Varying Mean (LVM)"  />     <Local_Mean_Property  value="Field_3:_thickness_or_lobe"  />     <Hard_Data  grid="Eros_data"   property="Thickness"  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="12" />     <Search_Ellipsoid  value="70 70 1  0 0 0" />    <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="1" filename =""   grid ="Smooth_sim_thick"  property ="Field_3:_thickness_or_lobe">  <LTI_type  function ="Power"  extreme ="-1.05"  omega ="3" />  <UTI_type  function ="No extrapolation"  extreme ="0"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="0" structures_count="1"  >    <structure_1  contribution="0.00200172"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>   
    end
else
    unconditional_sim=1
    if layer_flag==1 |layer_flag==0
        fprintf(fid2,'%s %s\n','RunGeostatAlgorithm ',['dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="Sim_Grid"  />     <Property_Name  value="',sim_name,'" />     <Nb_Realizations  value="1" />     <Seed  value="',num2str(ceil(100000+rand*100000)),'" />     <Kriging_Type  value="Kriging with Localy Varying Mean (LVM)"  />     <Local_Mean_Property  value="Thickness"  />     <Hard_Data  grid="',hard_data_grid,'   property="Thickness"  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="10" />     <Search_Ellipsoid  value="70 70 1  0 0 0" />    <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="0" filename =""   grid ="Sim_Grid"  property ="Thickness">  <LTI_type  function ="Power"  extreme ="0"  omega ="3" />  <UTI_type  function ="Power"  extreme ="',num2str(max_use),'"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="1e-006" structures_count="1"  >    <structure_1  contribution="',num2str(sill),'"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>']);
        % SKmean fprintf(fid2,'%s %s\n','RunGeostatAlgorithm ',['dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="Sim_Grid"  />     <Property_Name  value="',sim_name,'" />     <Nb_Realizations  value="1" />     <Seed  value="',num2str(ceil(100000+rand*100000)),'" />     <Kriging_Type  value="Simple Kriging (SK)"  />     <SK_mean  value="0.0" />     <Hard_Data  grid=""   property=""  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="10" />     <Search_Ellipsoid  value="70 70 1  0 0 0" />    <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="0" filename =""   grid ="Sim_Grid"  property ="Thickness">  <LTI_type  function ="Power"  extreme ="0"  omega ="3" />  <UTI_type  function ="Power"  extreme ="',num2str(max_use),'"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="1e-006" structures_count="1"  >    <structure_1  contribution="',num2str(sill),'"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>']);
        %                        RunGeostatAlgorithm  dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="Sim_Grid"  />     <Property_Name  value="dssim_data_map_1" />     <Nb_Realizations  value="1" />     <Seed  value="14071789" />     <Kriging_Type  value="Simple Kriging (SK)"  />     <SK_mean  value="0.0" />     <Hard_Data  grid="hard_grid_dssim	"   property="Thickness"  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="12" />                                <Search_Ellipsoid  value="70 70 1  0 0 0" />    <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="1" filename =""   grid ="Sim_Grid"  property ="Thickness">  <LTI_type  function ="Power"  extreme ="0"  omega ="3" />  <UTI_type  function ="Power"  extreme ="3.6369"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="1e-006" structures_count="1"  >    <structure_1  contribution="0.05309"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>   

        %                           RunGeostatAlgorithm     dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="IngFgu"  />     <Property_Name  value="test_lvm_1_uncond" />     <Nb_Realizations  value="1" />     <Seed  value="1407171" />                <Kriging_Type  value="Kriging with Localy Varying Mean (LVM)"  />     <Local_Mean_Property  value="Field_3:_thickness_or_lobe"  />     <Hard_Data  grid=""   property=""  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="12" />     <Search_Ellipsoid  value="70 70 1  0 0 0" />    <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="1" filename =""   grid ="IngFgu"  property ="Field_3:_thickness_or_lobe">  <LTI_type  function ="No extrapolation"  extreme ="-1.05"  omega ="3" />  <UTI_type  function ="Power"  extreme ="0.6"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="0" structures_count="1"  >    <structure_1  contribution="0.00685645"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>
    elseif layer_flag==2
        fprintf(fid2,'%s %s\n','RunGeostatAlgorithm ',['dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="Sim_Grid"  />     <Property_Name  value="',sim_name,'" />     <Nb_Realizations  value="1" />     <Seed  value="',num2str(ceil(100000+rand*100000)),'" />     <Kriging_Type  value="Kriging with Localy Varying Mean (LVM)"  />     <Local_Mean_Property  value="Thickness"  />     <Hard_Data  grid="',hard_data_grid,'   property="Thickness"  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="10" />     <Search_Ellipsoid  value="70 70 1  0 0 0" />     <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="0" filename =""   grid ="Sim_Grid"  property ="Thickness">  <LTI_type  function ="Power"  extreme ="',num2str(min_use),'"  omega ="3" />  <UTI_type  function ="No extrapolation"  extreme ="0"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="1e-006" structures_count="1"  >    <structure_1  contribution="',num2str(sill),'"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>']);
        % SKmean fprintf(fid2,'%s %s\n','RunGeostatAlgorithm ',['dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="Sim_Grid"  />     <Property_Name  value="',sim_name,'" />     <Nb_Realizations  value="1" />     <Seed  value="',num2str(ceil(100000+rand*100000)),'" />     <Kriging_Type  value="Simple Kriging (SK)"  />     <SK_mean  value="0.0" />     <Hard_Data  grid=""   property=""  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="10" />     <Search_Ellipsoid  value="70 70 1  0 0 0" />     <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="0" filename =""   grid ="Sim_Grid"  property ="Thickness">  <LTI_type  function ="Power"  extreme ="',num2str(min_use),'"  omega ="3" />  <UTI_type  function ="No extrapolation"  extreme ="0"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="1e-006" structures_count="1"  >    <structure_1  contribution="',num2str(sill),'"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>']);
                         %parfile: RunGeostatAlgorithm  dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="Sim_Grid"  />     <Property_Name  value="dssim_data_map_1" />     <Nb_Realizations  value="1" />     <Seed  value="153086" />     <Kriging_Type  value="Kriging with Localy Varying Mean (LVM)"  />     <Local_Mean_Property  value="Thickness"  />     <Hard_Data  grid="hard_grid_dssim	"   property="Thickness"  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="10" />     <Search_Ellipsoid  value="70 70 1  0 0 0" />    <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="1" filename =""   grid ="Sim_Grid"  property ="Thickness">  <LTI_type  function ="No extrapolation"  extreme ="0"  omega ="3" />  <UTI_type  function ="Power"  extreme ="3.6369"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="0" structures_count="1"  >    <structure_1  contribution="0.05309"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>
                         % SGEMS: RunGeostatAlgorithm  dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="Smooth_sim_thick"  />     <Property_Name  value="test_lvm_1" />     <Nb_Realizations  value="1" />     <Seed  value="14071791" />     <Kriging_Type  value="Kriging with Localy Varying Mean (LVM)"  />     <Local_Mean_Property  value="Field_3:_thickness_or_lobe"  />     <Hard_Data  grid="Hard_data_sst"   property="Thickness"  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="12" />     <Search_Ellipsoid  value="70 70 1  0 0 0" />    <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="1"   property ="Field_3:_thickness_or_lobe">  <LTI_type  function ="No extrapolation"  extreme ="0.001"  omega ="3" />  <UTI_type  function ="Power"  extreme ="3.5"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="0" structures_count="1"  >    <structure_1  contribution="0.0730641"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>   

                         % SGEMS: RunGeostatAlgorithm  dssim::/GeostatParamUtils/XML::<parameters>  <algorithm name="dssim" />     <Grid_Name  value="Smooth_sim_thick"  />     <Property_Name  value="test_lvm_1" />     <Nb_Realizations  value="1" />     <Seed  value="14071791" />     <Kriging_Type  value="Kriging with Localy Varying Mean (LVM)"  />     <Local_Mean_Property  value="Field_3:_thickness_or_lobe"  />     <Hard_Data  grid="Hard_data_sst"   property="Thickness"  />     <Assign_Hard_Data  value="1"  />     <Max_Conditioning_Data  value="12" />     <Search_Ellipsoid  value="70 70 1  0 0 0" />    <cdf_type  value="Soares"  />     <LN_mean  value="1" />     <LN_variance  value="1" />     <U_min  value="0" />     <U_max  value="1" />     <nonParamCdf  ref_on_file ="0"  ref_on_grid ="1"  break_ties ="1" filename =""   grid ="Smooth_sim_thick"  property ="Field_3:_thickness_or_lobe">  <LTI_type  function ="Power"  extreme ="0.001"  omega ="3" />  <UTI_type  function ="Power"  extreme ="3.5"  omega ="0.333" />  </nonParamCdf>    <is_local_correction  value="1"  />     <Variogram  nugget="0" structures_count="1"  >    <structure_1  contribution="0.0730641"  type="Gaussian"   >      <ranges max="500"  medium="500"  min="1"   />      <angles x="0"  y="0"  z="0"   />    </structure_1>  </Variogram>  </parameters>   

    end
end

fprintf(fid2,'%s %s\n','SaveGeostatGrid ',strcat(' Sim_Grid::',sgems_output_path_parfile,sim_filename,'::gslib::0::',strcat(sim_name,'__real0')));
fclose(fid2);


% Then run the algorithm
cd(sgems_code_path)
system('sgems dssim_parfile.log');
cd(matlab_code_path)

%read simulated grid
new_layerpath=strcat(sgems_output_path_normal,sim_filename);
dssim_thick_map=read_SGEMS_output(new_layerpath,g, 1);

% if cond==0  % if it's unconditional, had to simulate only the noise and now add it to the mean values
%     dssim_thick_map=dssim_thick_map+smooth_sim_thick;
% end

% figure;
% imagesc(g.x,g.y,dssim_thick_map);
% title('DSSIM of thicknesses');
% colorbar;
% xlabel('x [m]');
% ylabel('y [m]');
% axis xy
% axis equal
% hold on
% 
% if cond==1
%     for n=1:size(next_well_data,1)
%         if next_well_data(wnum,5)~=999  % if there's an interval left in the well
%             if next_well_data(n,4)==0   % fine-grained
%                 if next_well_data(n,7)==0|next_well_data(n,8)==-3   %bypass
%                     plot(g.dx*next_well_data(n,2),g.dy*next_well_data(n,3),'ko');
%                 else    % match (no bypass)
%                     plot(g.dx*next_well_data(n,2),g.dy*next_well_data(n,3),'wo');
%                 end
%             else    % if coarse-grained
%                 if next_well_data(n,7)==0|next_well_data(n,8)==-1   % bypass
%                     plot(g.dx*next_well_data(n,2),g.dy*next_well_data(n,3),'kd');
%                 else    % match (no bypass)
%                     plot(g.dx*next_well_data(n,2),g.dy*next_well_data(n,3),'yd');
%                 end
%             end
%         end
%     end
% end

