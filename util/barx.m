function barx(g)
% BARX create a stacked bar graph with negative data values on the negative 
ax1 = subplot(2,1,1,'XTickLabel',[]);
bar(g.*(g>0),'stacked')
axis tight
ax2 = subplot(2,1,2 );
bar(g.*(g<0),'stacked')
axis tight
lim1 = get(ax1,'YLim');
lim2 = get(ax2,'YLim');
pos = get(ax2,'position');
maxh = 1-2*pos(2);
posh = maxh*sum(abs(lim2))/sum(abs(lim1)+abs(lim2));
set(ax2,'position',[pos(1:3) posh])
set(ax1,'position',[pos(1) pos(2)+posh pos(3) maxh-posh])