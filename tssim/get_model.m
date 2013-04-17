% computes grid model for geoeas software (SGeMS) usage 
%[model_facies model_sec_var model_facies_geoeas model_sec_var_geoeas] = get_model(surface_cube,tracklobe,ny,nz,write_file,file_name,c,max_lobe_thickness_allowed,h)
% inputs:
% nz: number of blocks in z direction for discretization
% grid_block: grib block size (see get_model for computation) 
% write_file: flag that indicates if the file is written (filename)
%            the file is save in the Outputs folder
% c:parameter for distrubution of secondary variable ( 0 <=c<=1),ref.c=0.3
% max_lobe_thickness_allowed: user defined parameter (see in
%                             tssim)
% h: weight for local and max lobe thickness in the simulation( 0 <=c<=1)
%    ref. value = 0.3
% outputs:
%     model: Matlab 3D model
%     model_geoeas: 3D model file for geoeas software (SGeMS) usage
% no data = 1.2

% by Alejandro Leiva

function [model_facies model_sec_var model_facies_geoeas model_sec_var_geoeas] = get_model(surface_cube,tracklobe,ny,nz,write_file,file_name,c,max_lobe_thickness_allowed,h)
dirname = get_dirname;
[m n l]=size(surface_cube);
plots = find(tracklobe==-2);
if numel(plots)==0
    plots=l;
else
    plots = plots(1)-1;
end
d_max_h = max(max(surface_cube(:,:,plots)));
grid_block = d_max_h/nz;
model_facies = zeros(nz,n,ny);
model_sec_var = zeros(nz,n,ny);
for i=1:ny
    [model_facies(:,:,i) A zone_thickness] = get_model_slice(surface_cube,tracklobe,nz,i,grid_block,false,0);
    lower = find(A > 0 & A <= c);
    upper = find(A > c & A <= 1);
    A(lower) = A(lower)/c;
    A(upper) = (1-A(upper))/(1-c);
    if max_lobe_thickness_allowed~=0
       mid_range =  find(A > 0 & A <= 1);
       A(mid_range) = A(mid_range)*(1-h) + zone_thickness(mid_range)/max_lobe_thickness_allowed*h;
    end
    model_sec_var(:,:,i) =A;
        
    imagesc(model_facies(:,:,i));
    
    %imagesc(model_sec_var(:,:,i));
    set(gcf,'Color',[1 1 1]);scrsz = get(0,'ScreenSize');
    set(gcf,'Position',[100 scrsz(4)/3 scrsz(3)/1.2 scrsz(4)/2]);
    set(gca,'LineStyle','--','PlotBoxAspectRatio',[3.5 1 1]);colorbar;pause(0.001);
end
model_facies =flipdim(permute(model_facies,[3 2 1]),3);
model_sec_var =flipdim(permute(model_sec_var,[3 2 1]),3);
model_facies_geoeas = matlab2geoeas(model_facies);
model_sec_var_geoeas = matlab2geoeas(model_sec_var);
if write_file == true
    Input_path = strcat(dirname,'\code\Outputs\');
    fid =  fopen(strcat(Input_path,file_name),'w');
    fprintf(fid,'data\r\n');
    fprintf(fid,'%.0f\r\n',2);
    fprintf(fid,'facies\r\n');
    fprintf(fid,'sec_data\r\n');
    disp('Writing the geoeas file...');
    value = 1;
    for i = 1:length(model_facies_geoeas)
        if floor(10*i/length(model_facies_geoeas)) == value
            disp(strcat(num2str(100*i/length(model_facies_geoeas)),'% done'))
            value = value + 1;
        end
        fprintf(fid,'%6.2f %6.2f\r\n',model_facies_geoeas(i),model_sec_var_geoeas(i));
    end
    fclose(fid)
end