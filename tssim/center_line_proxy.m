% Considers channel center_line (A) and trend map from sediment source,
% matrix called B in the code.

function pbb_map = center_line_proxy(A)

binary = (A>0);
binary = bwdist(binary);
binary = 1-(1/max(max(binary))).*binary;
pbb_map = binary;


% result = result.*(A>0);
% x = find(A>0);
% j=min(result(x));
% result = (result-j).*(result>0/(1-j));





% B = trend_map_maker(g,center_point,b);
% [B_non_zeros_i B_non_zeros_j] = ind2sub(size(A),find(B~=0));
% pbb_map = zeros(size(A));
% new_A = A.*B;
% [non_zeros_i non_zeros_j] = ind2sub(size(A),find(new_A~=0));
% min_dist = zeros(size(B_non_zeros_i));
% for i = 1:size(B_non_zeros_i,1)
%     dist = sqrt((B_non_zeros_i(i)-non_zeros_i).^2+(B_non_zeros_j(i)-non_zeros_j).^2);
%     minim = min(dist);
%     min_dist(i) = minim;
% %     max2 = max(dist);
% %     if max2 > maxim
% %        maxim = max2; 
% %     end
%     pbb_map(B_non_zeros_i(i),B_non_zeros_j(i)) = minim; 
% end
% maxim = max(min_dist);
% pbb_map = 1-pbb_map/maxim;
% pbb_map(B == 0) = 0;
% pbb_map = pbb_map/sum(sum(pbb_map));
% imagesc(B);
% figure;imagesc(pbb_map);


%code that computes over the whole A - space (just for ppt)
% B = trend_map_maker(g,center_point,b)
% pbb_map = zeros(size(A));
% [zeros_i zeros_j] = ind2sub(size(A),find(A==0));
% [non_zeros_i non_zeros_j] = ind2sub(size(A),find(A~=0));
% 
% for i = 1:size(zeros_i,1)
%     dist = sqrt((zeros_i(i)-non_zeros_i).^2+(zeros_j(i)-non_zeros_j).^2);
%     dist = min(dist);
%     pbb_map(zeros_i(i),zeros_j(i)) = dist;
% end
% ppb_map = pbb_map.*B;
% pbb_map = pbb_map/sum(sum(pbb_map));
% imagesc(B);
% figure;imagesc(pbb_map);


