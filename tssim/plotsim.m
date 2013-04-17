function plotsim(surface_cube)
%Plots surface cube output from tssim
%function plotsim(surface_cube)

[m n l]=size(surface_cube);
plots = find(min(min(surface_cube))==0);
if numel(plots)==0
    plots=l;
end
rows = ceil(plots(1)/3);
for i=1:plots
    subplot(rows,3,i);imagesc(surface_cube(:,:,i));
end