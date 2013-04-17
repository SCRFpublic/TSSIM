%gives animation of the model slice
%code: 1: X-Z section 
%      2: Y-Z section
%n:number of slice in the particular direction
%     if code = 1, n<=size(surface,1)
%     if code = 2, n<=size(surface,2)
%delay:pause after every image shows up (ref: 0.1)
%by Alejandro Leiva
function model_movie(surface_cube,tracklobe,n,code,delay,hold_image)
close all;
if code==1 
    for i=1:n
        plot_model_slice(surface_cube,tracklobe,i,1,false);pause(delay);clc;
        if hold_image== false && i<n plot(1,1); end
    end
else
    for i=1:n
        plot_model_slice(surface_cube,tracklobe,i,2,false);pause(delay);clc;
        if hold_image== false && i<n plot(1,1); end
    end
end
