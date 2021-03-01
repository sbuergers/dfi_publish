function addtext( writethis )
% addtext adds the specified text in the top left corner of the current
% plot
xl    = xlim;
xdist = xl(2) - xl(1);
yl    = ylim;
ydist = yl(2) - yl(1);
text(xl(1)+0.1*xdist,yl(2)-0.1*ydist,writethis, 'color', [.7 .7 .7], 'fontsize', 8)

end

