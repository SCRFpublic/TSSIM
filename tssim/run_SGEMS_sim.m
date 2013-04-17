function [sim_layer]=run_SGEMS_sim(nTI,ti_name_prefix,lobe_num,cat_probs,matlab_output_path,sgems_output_path_parfile,sgems_output_path_normal,matlab_code_path,sgems_code_path,g2s_code_path,sim_name_root,hard_data_filename,use_pfields,pfield_filenames,servosystem,templateX,templateY,template_angle,tau,g,verbose);
% This function runs SGEMS in command mode from MATLAB
% [sim_layer]=run_SGEMS_sim(nTI,ti_name_prefix,lobe_num,cat_probs,matlab_output_path,sgems_output_path_parfile,sgems_output_path_normal,matlab_code_path,sgems_code_path,g2s_code_path,sim_name_root,hard_data_filename,use_pfields,pfield_filenames,servosystem,templateX,templateY,template_angle,tau,g,verbose);
% This function runs SGEMS in command mode and inputs grids by
% first converting them to Binary. Everything called to the command line
% through MATLAB

% Input:
%  nTI: number of TIs used
%  ti_name_prefix: prefix before the lobe number on the TIs to use
%  nlobe: number of lobes in each stack
%  lobe_num: the number of the current lobe to be simulated
%  cat_probs: categorical probabilities from the set of TI's: TI marginal
%  ti_file_path: path to the folder where the TI's were saved
%  sgems_output_path_parfile: path to the folder where the snesim simulation will be saved, with backwards slashes from normal
%  sgems_output_path_normal: same path as above, but with regular slashes
%  matlab_code_path: path to the matlab code folder
%  sgems_code_path: path to the sgems folder
%  g2s_code_path: path to geoeas2sgems - gslib to binary format conversion
%  sim_name_root: the name of the simulation. To this will be added the lobe_num and '__real0'
%  hard_data_filename: binary filename for hard data
%  use_pfields: 1 if using soft data, 0 if no
%  pfield_filenames: vector of binary filenames for the pfield for each facies
%  servosystem: parameter between 0 and 1 that controls how well the histogram is matched
%  templateX: size in cells of the snesim template in the x-direction
%  templateY: size in cells of the snesim template in the y-direction
%  tau: tau parameters as a string, e.g., '1 1'
%  g: simulation grid
%  verbose: whether or not to write what it's doing
  

% Written by: Holly Michael
% January 2008


% Here convert input pfield to binary and save it to the file name that is
% used in SGEMS.
% somehow 5 multigrids seems not to work with the pfields!
if use_pfields==1
    n_multigrids=3;
else
    n_multigrids=5;
end

% Extra check that cat_probs adds to 1!
%cat_probs
string_catprobs=num2str(cat_probs,6);
back_catprobs=str2num(string_catprobs);
cat_probs=back_catprobs;
cat_probs(end)=1.0-(sum(cat_probs)-cat_probs(end))

num_cat=length(cat_probs);  % number of facies categories
parfile_name='multiTI_parfile.log';
sim_name=strcat(sim_name_root,'_',num2str(lobe_num));
sim_filename=strcat(sim_name,'.dat');

sum_facies=g.nx*g.ny;

while (sum_facies==g.nx*g.ny);
    % Write the parameter log file:
    fid=fopen(strcat(sgems_code_path,parfile_name),'wt');

    if(fid==-1),
        fprintf('write_xydat_grid: Error opening file %s\n',parfile_name);
        error('cannot open file');
    end

    fprintf(fid,'%s%s%d%s%d%s\n','NewCartesianGrid ','TI_Grid::',g.nx,'::',g.ny,'::1::1::1::1.0::0::0::0');
    fprintf(fid,'%s%s%d%s%d%s\n','NewCartesianGrid ','Sim_Grid::',g.nx,'::',g.ny,'::1::1::1::1.0::0::0::0');

   % loop over each TI, each time change input and output file name, keep
    % track, call the set of TI's ti_name_set
    ti_name_set='';
    gridcount=0;

        for j=1:nTI
            ti_name=strcat(ti_name_prefix,num2str(lobe_num),'_',num2str(j));
            tifilename=strcat(ti_name,'.out');
            fprintf(fid,'%s %s\n','LoadObjectFromFile ',strcat(matlab_output_path,tifilename,'::All'));
            if gridcount==0
                fprintf(fid,'%s %s\n','CopyProperty ',strcat('grid::Field_3:_thickness_or_lobe::TI_Grid::',ti_name,'::0::0'));
            else
                fprintf(fid,'%s %s\n','CopyProperty ',strcat('grid_',num2str(gridcount),'::Field_3:_thickness_or_lobe::TI_Grid::',ti_name,'::0::0'));
            end
            gridcount=gridcount+1;
            %fprintf(fid,'%s %s\n','DeleteObjects ','grid');
            if (j==nTI)
                ti_name_set=strcat(ti_name_set,ti_name);
            else
                ti_name_set=strcat(ti_name_set,ti_name,';');
            end
        end
        
if use_pfields==1
    %have pfield filenames as a list: pfield_filenames
    pfield_filenames_str=cellstr(pfield_filenames);
    pfield_list='';
    for k=1:num_cat
        pfieldname=strcat('p',num2str(k-1));
        pfield_list=strcat(pfield_list,pfieldname,';');
        fprintf(fid,'%s %s\n','LoadObjectFromFile ',strcat(matlab_output_path,char(pfield_filenames_str(k)),'::All'));
        fprintf(fid,'%s %s\n','CopyProperty ',strcat('grid_',num2str(gridcount),'::Field_3:_thickness_or_lobe::Sim_Grid::',pfieldname,'::0::0'));
        gridcount=gridcount+1;
        %fprintf(fid,'%s %s\n','DeleteObjects ','grid');
    end
    fprintf(fid,'%s %s\n','LoadObjectFromFile ',strcat(matlab_output_path,hard_data_filename,'::All'));
    pfield_list;
else
    fprintf(fid,'%s %s\n','LoadObjectFromFile ',strcat(matlab_output_path,hard_data_filename,'::All'));
    pfield_list='none'
end

if templateX>templateY
    template_angle=template_angle+90;
    maxtemplate=templateX;
    mintemplate=templateY;
else
    maxtemplate=templateY;
    mintemplate=templateX;
end
    % then run the algorithm once.
    fprintf(fid,'%s %s\n','RunGeostatAlgorithm ',['Snesim_mti::/GeostatParamUtils/XML::<parameters>  <algorithm name="Snesim_mti" />     <GridSelector_Sim  value="Sim_Grid"  />     <Property_Name_Sim  value="',sim_name,'" />     <Nb_Realizations  value="1" />     <Seed  value="',num2str(ceil(100000+rand*100000)),'" />     <_ti  value="TI_Grid"  />     <_ti_list count="',num2str(nTI),'"   value="',ti_name_set,'"  />     <Nb_Facies  value="',num2str(num_cat,'%5.4d'),'" />     <Marginal_Cdf  value="',num2str(cat_probs,6),'" />     <Max_Cond  value="120" />     <Search_Ellipsoid  value="',num2str(maxtemplate),' ',num2str(mintemplate),' 1 ', num2str(template_angle), ' 0 0" />    <Hard_Data  grid="hard_grid"   property="Facies number"  />     <Use_ProbField  value="',num2str(use_pfields),'"  />     <ProbField_properties count="',num2str(num_cat),'"   value="',pfield_list,'"  />     <TauModelObject  value="',tau,'" />     <VerticalPropObject  value=""  />     <VerticalProperties count="0"   value=""  />     <Use_Affinity  value="0"  />     <Use_Rotation  value="0"  />     <Cmin  value="1" />     <Constraint_Marginal_ADVANCED  value="',num2str(servosystem),'" />     <resimulation_criterion  value="-1" />     <resimulation_iteration_nb  value="1" />     <Nb_Multigrids_ADVANCED  value="',num2str(n_multigrids),'" />     <Debug_Level  value="0" />     <Subgrid_choice  value="0"  />     <expand_isotropic  value="1"  />     <expand_anisotropic  value="0"  />     <aniso_factor  value="    " />     <Region_Indicator_Prop  value=""  />     <Active_Region_Code  value="" />     <Use_Previous_Simulation  value="0"  />     <Use_Region  value="0"  />   </parameters>']);
    fprintf(fid,'%s %s\n','SaveGeostatGrid ',strcat(' Sim_Grid::',sgems_output_path_parfile,sim_filename,'::gslib::0::',strcat(sim_name,'__real0')));

    fclose(fid);
    % Then run the algorithm
    cd(sgems_code_path)
    system('sgems multiTI_parfile.log')
    cd(matlab_code_path)

    %read simulated grid
    new_layerpath=strcat(sgems_output_path_normal,sim_filename);
    sim_layer=read_SGEMS_output(new_layerpath,g, verbose);

    % Here put a test: if sum(new_layer)=g.nx*g.ny, then run sgems again.

    sum_facies=sum(sum(sim_layer));  % If the simulation has no lobe, everything will be 1, and the sum will equal the number of grid cells

end % This ends the while loop that checks to see that the simulation has a lobe in it.

% sim_marginal=zeros(1,num_cat);
% for m=1:num_cat
%     sim_marg_tot(m)=length(find(sim_layer==(m-1)));
% end
% sim_marginal=sim_marg_tot/g.nx/g.ny
% input_marginal=cat_probs    


% figure;
% imagesc(g.x,g.y,sim_layer);
% title('Snesim-simulated channel-lobe thickness map');
% colorbar;
% xlabel('x [m]');
% ylabel('y [m]');
% axis xy
% axis equal


