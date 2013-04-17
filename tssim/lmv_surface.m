% smoothes out surface E using local moving windows of windows_size.
% inputs:
%     E: DEM.
%     Windows_size: Windows size used for lmv.
%     g: grid dimension.
% outputs: 
%     lmved_surface: smoothened surfaced.
% by Alejandro Leiva
    
function lmved_surface = lmv_surface(E,windows_size,g)
E1=zeros(g.ly/g.dy,g.lx/g.dx,2);
E1(:,:,1)=E;
E2=smooth3(E1,'box',[windows_size windows_size 1]);
lmved_surface=E2(:,:,1);