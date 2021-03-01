function fh = plot_correlation(mat1, mat2, xlbl, ylbl)
    fh = figure('color', 'w', 'position', [50 50 900, 250]);
    cond_labels = {'0S', '1S', '2S'};
    for icond = 1:3
        subplot(1,3,icond)
        [spearRho, pval] = corr(mat1(:,icond), mat2(:,icond), 'type', 'Spearman', 'rows', 'complete');
        [~,b1,b0] = regression(mat1(:,icond), mat2(:,icond), 'one');
        plot(mat1(:,icond), mat2(:,icond), 'bo', 'markersize', 10); hold on
        l1 = min([mat1(:,icond); mat2(:,icond)]); 
        l2 = max([mat1(:,icond); mat2(:,icond)]);
        line([l1  l2], [b0+b1*l1 b0+b1*l2], 'color', 'm');
        xlim([l1, l2]); ylim([l1, l2]);
        addtext(sprintf('\nrho = %.2f\np-value = %.5f\nN = %.2f', spearRho, pval, 20));
        xlabel(xlbl);
        ylabel(ylbl);
        title(cond_labels{icond});
    end
end