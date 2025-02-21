function [] = plot_rt_probs(P_matrix)
    I = mat2gray((P_matrix.^0.175)');
    I = 1-I;
    %I = flipud(I);
    figure();
    imshow(I);
    set(gca, 'YDir','normal');
    ylim([0 500])
    axis on;
    colorbar('Ticks',[0 1],'Ticklabels',[1 0],'Direction','reverse');
    title('Run length Densities'); xlabel('Time index');ylabel('Run Length');
    figure;
    [~, ind] = max(P_matrix,[],2);
    plot(ind);
    title('Predicted run length');xlabel('Time index');ylabel('Run Length');
end