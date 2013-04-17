%gets grid model for slice number nslice
%nz: number of blocks in z direction for discretization
%grid_block: grid block size (see get_model for computation) 
%log_or_sec: plotting --> 1(log), 0(secondary data)
%[log_facies log_sec_data zone_thickness] = get_model_slice(surface_cube,tracklobe,nz,nslice,grid_block,plot_figure,log_or_sec)

%by Alejandro Leiva
function [log_facies log_sec_data zone_thickness] = get_model_slice(surface_cube,tracklobe,nz,nslice,grid_block,plot_figure,log_or_sec)
[m n l]=size(surface_cube);
slice = reshape(surface_cube(nslice,:,:),n,l)';
log_facies = zeros(nz,n);
log_sec_data = zeros(nz,n);
zone_thickness = zeros(nz,n);
for j=1:n
    altitude = grid_block;
    for i=1:nz
        pos = find(slice(:,j)>altitude);
        if numel(pos)==0 %air
            log_facies(i:nz,j) = 1.2;
            log_sec_data(i:nz,j) = 1.2;
            zone_thickness(i:nz,j) = 1.2;
            break;
        end
        facie = tracklobe(pos(1));
        if facie == 0%shale
            element = 0;
            log_sec_data(i,j) = 0;
        end
        if facie > 0%sand
            element = 3;
            
            if pos(1)-1 == 0
                dh = altitude - 0;
                dz = slice(pos(1),j)-0;
            else 
                dh = altitude - slice(pos(1)-1,j);
                dz = slice(pos(1),j)-slice(pos(1)-1,j);
            end
            log_sec_data(i,j) = dh/dz;
            zone_thickness(i,j) = dz;
        end
        log_facies(i,j) = element;
        
        if slice(1,j)>altitude%underneath surface
            log_facies(i,j) = 1.2;
            log_sec_data(i,j) = 1.2;
            zone_thickness(i,j) = 1.2;
        end
        altitude = altitude + grid_block;      
    end
end
log_facies = flipud(log_facies);
log_sec_data = flipud(log_sec_data);
zone_thickness = flipud(zone_thickness);
if plot_figure == true
    if log_or_sec == 0
        figure;imagesc(log_sec_data);colorbar;
    else
        figure;imagesc(log_facies);colorbar;
    end
end
