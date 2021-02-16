close all;
P_matrix = P_runlength;
% mat = log(ones(size(P_matrix))-P_matrix);
% for k = 1:length(mat)
%     mat(k,:) = mat(k,:)/sum(mat(k,:));
% end
I = (P_matrix.^0.2)';
I = 1-I;
I = flipud(I);
figure();

imshow(I);
axis on;
colorbar('Ticks',[0 1],'Ticklabels',[1 0],'Direction','reverse');
title('Run length Densities'); xlabel('Time index');ylabel('Run Length');
figure;
[~, ind] = max(P_matrix,[],2);
plot(ind);
title('Predicted run length');xlabel('Time index');ylabel('Run Length');