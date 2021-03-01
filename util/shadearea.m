function shadearea(xrange)
%Highlight area between two x values in plot
% xrange specifies those x boundaries
yl = ylim;
harea = area(xrange, [yl(1), yl(1)], yl(2), 'LineStyle', 'none');
set(harea, 'FaceColor', 'b')
alpha(0.25);

end

