function plotspecs
    % general specifications you want to apply quickly to a plot
    % do not change legend objects...
    %legh = findobj(gcf,'Type','axes','Tag','legend');
    %leg_old_fs = get(legh, 'FontSize');
    set(findall(gcf,'-property','FontSize'),'FontSize',12)
    set(findall(gcf,'-property','LineWidth'), 'LineWidth', 1.5)
    set(gcf,'color','w');
    set(findall(gcf,'type','axes'),'YDir','normal');
    %set(legh, 'FontSize', leg_old_fs);
end
