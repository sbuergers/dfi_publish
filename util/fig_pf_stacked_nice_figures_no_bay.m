function fh = fig_pf_stacked_nice_figures_no_bay(pffit, ylab, ttl, distr, position)
% Plot psychometric functions of multiple participants next to each other.
% Conditions are plotted below each other.
%
% INPUT: s      =  structure from dfi paradigm
%        pffit  =  output from dfi_fit_pf in the form {participant}{condition}
%        ylab   =  cell array with labels for conditions
%        ttl    =  column titles in cell array
%        supttl =  figure title
%         distr =  which psychometric function type to use (e.g. 'Normal')
%      position =  position and dimensions of the figure
%
% OUTPUT: fh    =  the figure handle
%


if exist('distr', 'var')
    switch distr
        case 'Weibull', PF = @PAL_Weibull;
        case 'Normal',  PF = @PAL_CumulativeNormal;
        case 'Gumbel',  PF = @PAL_Gumbel;           % log Weibull
        case 'Quick',   PF = @PAL_Quick;
        case 'logQuick',PF = @PAL_logQuick;
        case 'HypSec',  PF = @PAL_HyperbolicSecant;
        case 'logist',  PF = @PAL_Logistic;    
    end
else
    distr = 'Normal';
    PF = @PAL_CumulativeNormal;
end 


% Plot previous figures next to each other
% size
if ~exist('position', 'var')
    positionXY = [0, 0, 600*numel(pffit), 1020];
else
    positionXY = position;
end



ni = numel(pffit{1}); 
nj = numel(pffit);
fh = figure('color', [1 1 1], 'Position', positionXY);
base.left   = 100;% pixels from left
base.bottom = 80; % pixels from bottom
jump.left   = 10; % pixels between plots laterally
jump.bottom = 10; % pixels between plots vertically
figsize.x   = (positionXY(3)-base.left*2  -jump.left*(nj-1))/nj;
figsize.y   = (positionXY(4)-base.bottom*2-jump.bottom*(ni-1))/ni;

col.points = 'k';
col.funML  = [0.3   0.8   0.3];
col.ipML   = [0.1   0.6   0.1];
col.funBA  = [0.8   0.3   0.3];
col.ipBA   = [0.6 0.1 0.1];

lw   = 4;  % PF line widths
dtsz = 10; % size of PF dots

for i = 1:ni % condition
    for j = 1:nj % subj
%         try
            
            soa  = pffit{j}{i}.soa;
            soal = round((soa*1000))/1000; % round for plotting to two decimal places

            % determine plot id and position
            pid = (i-1)*nj+j;
            subplot(ni,nj,pid, 'Position', [(base.left+(j-1)*figsize.x+(j-1)*jump.left)/positionXY(3), ...
                            (positionXY(4)-base.bottom-(i-1)*jump.bottom-(i)*figsize.y)/positionXY(4), ...
                             figsize.x/positionXY(3), figsize.y/positionXY(4)]);
            % Plot data-points
            plot(soa, pffit{j}{i}.perCor, 'ko', 'MarkerSize', dtsz, 'MarkerFaceColor', col.points); 
            % Plot PF
            hold on
            soaHR = min(soa):max(soa)/1000:max(soa);
            PFfitML = PF(pffit{j}{i}.par,    soaHR);
            hML = plot(soaHR, PFfitML, '--', 'Color', col.funML, 'LineWidth', lw); %'g-'
            % plot settings
            set(gca, 'Xtick', soa); yl = [-0.05, 1.05]; ylim(yl);
            set(gca, 'Xlim', [0 max(soa)+min(soa)])
            line([pffit{j}{i}.par(1) pffit{j}{i}.par(1)], [yl], 'color', col.ipML, 'LineWidth', 1.5) % threshold / IP ML
            % only set for first row
            if i == 1, title(ttl{j}); end
            % only set for last row
            if i == ni
                xlabel('SOA'); 
                set(gca, 'XtickLabel', round(soa*1000)/1000); rotateXLabels(gca,75); 
            end
            % set for first two rows
            if i ~= ni, set(gca, 'XTickLabel', []); end
            % set for first column
            if j == 1, ylabel(sprintf([ylab{i}, '\nAccuracy'])); end
            % set for all columns but the first
            if j ~= 1, set(gca, 'YTickLabel', []); end
            
            gridxy(get(gca,'xtick'),get(gca,'ytick'),'color',[.6 .6 .6],'linewidth',1, 'linestyle', ':')
            if i == 1
                set(gca, 'color', [0.9  1 0.9 ])
            elseif i == 2
                set(gca, 'color', [0.9  0.9  1])
            elseif i == 3 
                set(gca, 'color', [1 0.9  0.9 ])
            end
            
            % add parameter values as text
            text(max(soa)+0.02, 0.24, sprintf('treshold (%.2f)\nslope (%.2f)\nguess-rate (%.2f)\nlapse-rate (%.2f)\neta (%.2f)', ...
                                              pffit{j}{i}.par(1), ...
                                              pffit{j}{i}.par(2), ...
                                              pffit{j}{i}.par(3), ...
                                              pffit{j}{i}.par(4), ...
                                              pffit{j}{i}.eta), ...
                                              'FontSize', 10, 'HorizontalAlignment', ...
                                              'right', 'VerticalAlignment' , 'middle');
%         catch
%             fprintf('\n\nNot able to plot this particular condition...\n');
%         end
    end % j
end % i



% eof 

















