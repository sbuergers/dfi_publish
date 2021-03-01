% Create figures of pre-stimulus frequency within subject
% analysis for both dprime and bias, and their bayes factors.
%
% Parent script(s): 
%   ?
%   
% Children script(s): 
%   None
%
% Sibling script(s):
%   None
%
%
% DETAILS
%
% None
% 
% ---
% Steffen Buergers, sbuergers@gmail.com,
% Last modified Feb. 2021



%% Setup

clear all
close all
clc

try
    addpath(genpath('dfi'))
catch
    warning('Cannot find dfi folder')
end
dfi_startup

data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
fig_dir  = fullfile('dfi_experiment_figures');

try
    addpath(fullfile('fieldtrip-20160816'))
catch
    warning('Cannot find fieldtrip folder')
end

% useful variables
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N = length(subjvect);
condvect = {'v2', 'Fus', 'Fis', 'av2'};
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]};  % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]};  % see2

% figure rendering
useopengl = false;



%% Figure 1 (dprime)

fh1 = figure('color', [1 1 1], 'Position', [0, 0, 427, 705]);
ha = tight_subplot(6, 4,[0.02 0.02],[0.02],[0.02]);

%% Row 1 (yn_intermsoas, sensor)
clear fslide_GA tif sem_f1 sem_f2 stat

xl = [-0.6 -0.1];
yl = [0.9, 2.4];

ids = 1:3;
for icond = 1:3
    axes(ha(ids(icond)));
    
    % Load in figure data
    load(fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\fslide\yesno\dprime\no_regress', ...
            sprintf('figure_data_PO4_O2_PO8_%s.mat', condvect{icond})));
        
    % Which clusters are reliable?
    % Make a vector of all p-values associated with the clusters from ft_freqstatistics.
    if ~isfield(stat, 'posclusters')
        pos_cluster_pvals = [];
        pos = [];
    elseif isempty(stat.posclusters)
        pos_cluster_pvals = [];
        pos = [];
    else
        pos_cluster_pvals = [stat.posclusters(:).prob];
        % Then, find which clusters are significant, outputting their indices as held in stat.posclusters
        pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
        % make a boolean matrix of which (freq, time)-pairs are part of a significant cluster
        pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
    end
    
    % and now for the negative clusters...
    if ~isfield(stat, 'negclusters')
        neg_cluster_pvals = [];
        neg = [];
    elseif isempty(stat.negclusters)
        neg_cluster_pvals = [];
        neg = [];
    else
        neg_cluster_pvals = [stat.negclusters(:).prob];
        neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
        neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
    end
    
    % use pchip to interpolate data and make the figure less edgy
%     >> Y = [1,1,2,1,1];
%     >> X = [0,1,3,4,5];
%     >> Xi = 0:0.1:5;
%     >> Yi = pchip(X,Y,Xi);
%     >> plot(Xi,Yi,X,Y)
    
    time_bool = tif >= -0.602 & tif <= -0.1;
    fslide_GA = fslide_GA(time_bool,:);
    tif = stat.time(stat.time <= -0.1);
    tif_chip = tif(1):0.01:tif(end);
    shadedErrorBar(tif_chip, pchip(tif, fslide_GA(: ,2), tif_chip), pchip(tif, sem_f2_w', tif_chip), {'-b','markerfacecolor', col2{icond}}, useopengl); hold on
    shadedErrorBar(tif_chip, pchip(tif, fslide_GA(: ,1), tif_chip), pchip(tif, sem_f1_w', tif_chip), {'-r','markerfacecolor', col1{icond}}, useopengl);
    % accentuate actual GA time courses
    plot(tif_chip, pchip(tif, fslide_GA(: ,2), tif_chip),  'color', col2{icond}, 'linewidth', 1.5); hold on; % correct
    plot(tif_chip, pchip(tif, fslide_GA(: ,1), tif_chip),  'color', col1{icond}, 'linewidth', 1.5);          % incorrect
    % Include positive clusters in figure
    X = [tif(pos),fliplr(tif(pos))];
    Y = [fslide_GA(pos,2)', fliplr(fslide_GA(pos,1)')];
    fill(X,Y,[0.8, 0.8, 0.8], 'linewidth', 0.5);
    % Include negative clusters in figure
    X = [tif(neg),fliplr(tif(neg))];
    Y = [fslide_GA(neg,2)', fliplr(fslide_GA(neg,1)')];
    fill(X,Y,[0.8, 0.8, 0.8], 'linewidth', 0.5);
    % add grid lines
    %gridxy(get(gca, 'xtick'), get(gca, 'ytick'), 'color', [0.9, 0.9, 0.9])
    xlim(xl)
    xticks = -0.6:0.1:0;  set(gca, 'XTick', xticks);
    ylim(yl)
    yticks = 0:0.5:100;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    set(gca,'XColor','k','YColor','k')
    %set(gca, 'linewidth', 2)
    
end % condition loop

%% Row 4 (ynt sensor)
clear fslide_GA tif sem_f1 sem_f2 stat

xl = [-0.6 -0.1];
yl = [0.9, 2.4];

ids = 13:16;
for icond = 1:3
    axes(ha(ids(icond)));
    
    % Load in figure data
    load(fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\fslide\ynt\dPrime\no_regress', ...
            sprintf('figure_data_PO4_O2_PO8_%s.mat', condvect{icond})));
        
    % Which clusters are reliable?
    % Make a vector of all p-values associated with the clusters from ft_freqstatistics.
    if ~isfield(stat, 'posclusters')
        pos_cluster_pvals = [];
        pos = [];
    elseif isempty(stat.posclusters)
        pos_cluster_pvals = [];
        pos = [];
    else
        pos_cluster_pvals = [stat.posclusters(:).prob];
        % Then, find which clusters are significant, outputting their indices as held in stat.posclusters
        pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
        % make a boolean matrix of which (freq, time)-pairs are part of a significant cluster
        pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
    end
    
    % and now for the negative clusters...
    if ~isfield(stat, 'negclusters')
        neg_cluster_pvals = [];
        neg = [];
    elseif isempty(stat.negclusters)
        neg_cluster_pvals = [];
        neg = [];
    else
        neg_cluster_pvals = [stat.negclusters(:).prob];
        neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
        neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
    end
    
    % Plot results (gray area denotes significant cluster)
    %opengl('OpenGLDockingBug',1)
    %add error bars
    time_bool = tif >= -0.602 & tif <= -0.1;
    fslide_GA = fslide_GA(time_bool,:);
    tif = tif(time_bool);
    shadedErrorBar(tif_chip, pchip(tif, fslide_GA(: ,2), tif_chip), pchip(tif, sem_f2_w', tif_chip), {'-b','markerfacecolor', col2{icond}}, useopengl); hold on
    shadedErrorBar(tif_chip, pchip(tif, fslide_GA(: ,1), tif_chip), pchip(tif, sem_f1_w', tif_chip), {'-r','markerfacecolor', col1{icond}}, useopengl);
    % accentuate actual GA time courses
    plot(tif_chip, pchip(tif, fslide_GA(: ,2), tif_chip),  'color', col2{icond}, 'linewidth', 1.5); hold on; % correct
    plot(tif_chip, pchip(tif, fslide_GA(: ,1), tif_chip),  'color', col1{icond}, 'linewidth', 1.5);          % incorrect
    % Include positive clusters in figure
    X = [tif(pos),fliplr(tif(pos))];
    Y = [fslide_GA(pos,2)', fliplr(fslide_GA(pos,1)')];
    fill(X,Y,[0.8, 0.8, 0.8], 'linewidth', 0.5);
    % Include negative clusters in figure
    X = [tif(neg),fliplr(tif(neg))];
    Y = [fslide_GA(neg,2)', fliplr(fslide_GA(neg,1)')];
    fill(X,Y,[0.8, 0.8, 0.8], 'linewidth', 0.5);
    % add grid lines
    %gridxy(get(gca, 'xtick'), get(gca, 'ytick'), 'color', [0.9, 0.9, 0.9])
    xlim(xl)
    xticks = -0.6:0.1:0;  set(gca, 'XTick', xticks);
    ylim(yl)
    yticks = -100:0.5:100;;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    
end % condition loop


export_fig(fh1, 'D:\dfi_experiment_figures\Paper_figures\iAF\iAF_within\sd_params\fslide_dprime_no_regress', '-tiff', '-m2.5');
fh.Renderer = 'painters'; 
saveas(fh1,fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\iAF_within\sd_params\fslide_dprime_no_regress_svg.svg'))
close all


%% Figure 2 (Bayes factors - dprime)

fh2 = figure('color', [1 1 1], 'Position', [0, 0, 427, 400]);
ha = tight_subplot(6, 4, [0.02 0.02], [0.02], [0.02]);

%% Row 1 (yn_intermsoas, sensor)
clear fslide_GA tif sem_f1 sem_f2 stat

xl = [-0.6 -0.1];
yl = [-1.25 1.75];

ids = 1:3;
for icond = 1:3
    axes(ha(ids(icond)));

    load(fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\fslide\yesno\dprime\no_regress', ...
        sprintf('Bayes_factors_PO4_O2_PO8_%s.mat', condvect{icond})), 'bf');

    tif = linspace(-0.6016, -0.1211, length(bf));
    tif_chip = tif(1):0.01:tif(end);
    % add jeffrey's substantial evidence lines
    plot(tif, repmat(-0.5, size(tif)), 'color', [0.75 0.25 0.75]); hold on
    plot(tif, repmat(+0.5, size(tif)), 'color', [0.75 0.25 0.75]);
    % plot bayes factor time series
    plot(tif_chip, pchip(tif, log10(bf), tif_chip),  'color', 'k', 'linewidth', 1); 
    xlim(xl)
    xticks = -0.6:0.1:0;  set(gca, 'XTick', xticks);
    ylim(yl)
    yticks = -2:0.5:2;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    set(gca,'XColor','k','YColor','k')
    %set(gca, 'linewidth', 2)
    
end % condition loop

%% Row 4 (ynt sensor)
clear fslide_GA tif sem_f1 sem_f2 stat

xl = [-0.6 -0.1];

ids = 13:16;
for icond = 1:3
    axes(ha(ids(icond)));
    
    % Save data to remake figures later (for paper)
    load(fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\fslide\ynt\dPrime\no_regress', ...
        sprintf('Bayes_factors_PO4_O2_PO8_%s.mat', condvect{icond})), 'bf');
    
    tif = linspace(-0.6016, -0.1211, length(bf));
    tif_chip = tif(1):0.01:tif(end);
    % add jeffrey's substantial evidence lines
    plot(tif, repmat(-0.5, size(tif)), 'color', [0.75 0.25 0.75]); hold on
    plot(tif, repmat(+0.5, size(tif)), 'color', [0.75 0.25 0.75]);
    % plot bayes factor time series
    plot(tif_chip, pchip(tif, log10(bf), tif_chip),  'color', 'k', 'linewidth', 1); 
    xlim(xl)
    xticks = -0.6:0.1:0;  set(gca, 'XTick', xticks);
    ylim(yl)
    yticks = -2:0.5:2;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    
end % condition loop


export_fig(fh2, 'D:\dfi_experiment_figures\Paper_figures\iAF\iAF_within\sd_params\fslide_dprime_no_regress_bf', '-eps', '-tiff', '-m2.5');
plot2svg('D:\dfi_experiment_figures\Paper_figures\iAF\iAF_within\sd_params\fslide_dprime_no_regress_bf_svg',fh2)









%% Figure 3 (criterion)

fh1 = figure('color', [1 1 1], 'Position', [0, 0, 427, 705]);
ha = tight_subplot(6, 4,[0.02 0.02],[0.02],[0.02]);

%% Row 1 (yn_intermsoas, sensor)
clear fslide_GA tif sem_f1 sem_f2 stat

xl = [-0.6 -0.1];
yl = [-1.25, 1.25];

ids = 1:3;
for icond = 1:3
    axes(ha(ids(icond)));
    
    % Load in figure data
    load(fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\fslide\yesno\criterion\no_regress', ...
            sprintf('figure_data_PO4_O2_PO8_%s.mat', condvect{icond})));
        
    % Make sure bias is consistent between conditions (i.e. 2F is always
    % signal, and 1F is always noise). I computed this differently for the
    % double flash illusion, so now simply flip the sign to make it so. 
    if icond == 3
        fslide_GA = -fslide_GA;
    end
    
    % Which clusters are reliable?
    % Make a vector of all p-values associated with the clusters from ft_freqstatistics.
    if ~isfield(stat, 'posclusters')
        pos_cluster_pvals = [];
        pos = [];
    elseif isempty(stat.posclusters)
        pos_cluster_pvals = [];
        pos = [];
    else
        pos_cluster_pvals = [stat.posclusters(:).prob];
        % Then, find which clusters are significant, outputting their indices as held in stat.posclusters
        pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
        % make a boolean matrix of which (freq, time)-pairs are part of a significant cluster
        pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
    end
    
    % and now for the negative clusters...
    if ~isfield(stat, 'negclusters')
        neg_cluster_pvals = [];
        neg = [];
    elseif isempty(stat.negclusters)
        neg_cluster_pvals = [];
        neg = [];
    else
        neg_cluster_pvals = [stat.negclusters(:).prob];
        neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
        neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
    end
    
    % use pchip to interpolate data and make the figure less edgy
%     >> Y = [1,1,2,1,1];
%     >> X = [0,1,3,4,5];
%     >> Xi = 0:0.1:5;
%     >> Yi = pchip(X,Y,Xi);
%     >> plot(Xi,Yi,X,Y)
    
    time_bool = tif >= -0.602 & tif <= -0.1;
    fslide_GA = fslide_GA(time_bool,:);
    tif = stat.time(stat.time <= -0.1);
    tif_chip = tif(1):0.01:tif(end);
    shadedErrorBar(tif_chip, pchip(tif, fslide_GA(: ,2), tif_chip), pchip(tif, sem_f2_w', tif_chip), {'-b','markerfacecolor', col2{icond}}, useopengl); hold on
    shadedErrorBar(tif_chip, pchip(tif, fslide_GA(: ,1), tif_chip), pchip(tif, sem_f1_w', tif_chip), {'-r','markerfacecolor', col1{icond}}, useopengl);
    % accentuate actual GA time courses
    plot(tif_chip, pchip(tif, fslide_GA(: ,2), tif_chip),  'color', col2{icond}, 'linewidth', 1.5); hold on; % correct
    plot(tif_chip, pchip(tif, fslide_GA(: ,1), tif_chip),  'color', col1{icond}, 'linewidth', 1.5);          % incorrect
    % Include positive clusters in figure
    X = [tif(pos),fliplr(tif(pos))];
    Y = [fslide_GA(pos,2)', fliplr(fslide_GA(pos,1)')];
    fill(X,Y,[0.8, 0.8, 0.8], 'linewidth', 0.5);
    % Include negative clusters in figure
    X = [tif(neg),fliplr(tif(neg))];
    Y = [fslide_GA(neg,2)', fliplr(fslide_GA(neg,1)')];
    fill(X,Y,[0.8, 0.8, 0.8], 'linewidth', 0.5);
    % add grid lines
    %gridxy(get(gca, 'xtick'), get(gca, 'ytick'), 'color', [0.9, 0.9, 0.9])
    xlim(xl)
    xticks = -0.6:0.1:0;  set(gca, 'XTick', xticks);
    ylim(yl)
    yticks = -100:0.5:100;;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    set(gca,'XColor','k','YColor','k')
    %set(gca, 'linewidth', 2)
    
end % condition loop

%% Row 4 (ynt sensor)
clear fslide_GA tif sem_f1 sem_f2 stat

xl = [-0.6 -0.1];
yl = [-1.25, 1.25];

ids = 13:16;
for icond = 1:3
    axes(ha(ids(icond)));
    
    % Load in figure data
    load(fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\fslide\ynt\criterion\no_regress', ...
            sprintf('figure_data_PO4_O2_PO8_%s.mat', condvect{icond})));
    
    % Make sure bias is consistent between conditions (i.e. 2F is always
    % signal, and 1F is always noise). I computed this differently for the
    % double flash illusion, so now simply flip the sign to make it so. 
    if icond == 3
        fslide_GA = -fslide_GA;
    end
    
    % Which clusters are reliable?
    % Make a vector of all p-values associated with the clusters from ft_freqstatistics.
    if ~isfield(stat, 'posclusters')
        pos_cluster_pvals = [];
        pos = [];
    elseif isempty(stat.posclusters)
        pos_cluster_pvals = [];
        pos = [];
    else
        pos_cluster_pvals = [stat.posclusters(:).prob];
        % Then, find which clusters are significant, outputting their indices as held in stat.posclusters
        pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
        % make a boolean matrix of which (freq, time)-pairs are part of a significant cluster
        pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
    end
    
    % and now for the negative clusters...
    if ~isfield(stat, 'negclusters')
        neg_cluster_pvals = [];
        neg = [];
    elseif isempty(stat.negclusters)
        neg_cluster_pvals = [];
        neg = [];
    else
        neg_cluster_pvals = [stat.negclusters(:).prob];
        neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
        neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
    end
    
    % Plot results (gray area denotes significant cluster)
    %opengl('OpenGLDockingBug',1)
    %add error bars
    time_bool = tif >= -0.602 & tif <= -0.1;
    fslide_GA = fslide_GA(time_bool,:);
    tif = tif(time_bool);
    shadedErrorBar(tif_chip, pchip(tif, fslide_GA(: ,2), tif_chip), pchip(tif, sem_f2_w', tif_chip), {'-b','markerfacecolor', col2{icond}}, useopengl); hold on
    shadedErrorBar(tif_chip, pchip(tif, fslide_GA(: ,1), tif_chip), pchip(tif, sem_f1_w', tif_chip), {'-r','markerfacecolor', col1{icond}}, useopengl);
    % accentuate actual GA time courses
    plot(tif_chip, pchip(tif, fslide_GA(: ,2), tif_chip),  'color', col2{icond}, 'linewidth', 1.5); hold on; % correct
    plot(tif_chip, pchip(tif, fslide_GA(: ,1), tif_chip),  'color', col1{icond}, 'linewidth', 1.5);          % incorrect
    % Include positive clusters in figure
    X = [tif(pos),fliplr(tif(pos))];
    Y = [fslide_GA(pos,2)', fliplr(fslide_GA(pos,1)')];
    fill(X,Y,[0.8, 0.8, 0.8], 'linewidth', 0.5);
    % Include negative clusters in figure
    X = [tif(neg),fliplr(tif(neg))];
    Y = [fslide_GA(neg,2)', fliplr(fslide_GA(neg,1)')];
    fill(X,Y,[0.8, 0.8, 0.8], 'linewidth', 0.5);
    % add grid lines
    %gridxy(get(gca, 'xtick'), get(gca, 'ytick'), 'color', [0.9, 0.9, 0.9])
    xlim(xl)
    xticks = -0.6:0.1:0;  set(gca, 'XTick', xticks);
    ylim(yl)
    yticks = -100:0.5:100;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    
end % condition loop


export_fig(fh1, 'D:\dfi_experiment_figures\Paper_figures\iAF\iAF_within\sd_params\fslide_criterion_no_regress_dfi_bug_fixed', '-tiff', '-m2.5');
fh.Renderer = 'painters'; 
saveas(fh1,fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\iAF_within\sd_params\fslide_criterion_no_regress_dfi_bug_fixed.svg'))
close all


%% Figure 4 (Bayes factors - criterion)

fh2 = figure('color', [1 1 1], 'Position', [0, 0, 427, 400]);
ha = tight_subplot(6, 4,[0.02 0.02],[0.02],[0.02]);

%% Row 1 (yn_intermsoas, sensor)
clear fslide_GA tif sem_f1 sem_f2 stat

xl = [-0.6 -0.1];
yl = [-1.25 1.75];

ids = 1:3;
for icond = 1:3
    axes(ha(ids(icond)));

    load(fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\fslide\yesno\criterion\no_regress', ...
        sprintf('Bayes_factors_PO4_O2_PO8_%s.mat', condvect{icond})), 'bf');

    tif = linspace(-0.6016, -0.1211, length(bf));
    tif_chip = tif(1):0.01:tif(end);
    % add jeffrey's substantial evidence lines
    plot(tif, repmat(-0.5, size(tif)), 'color', [0.75 0.25 0.75]); hold on
    plot(tif, repmat(+0.5, size(tif)), 'color', [0.75 0.25 0.75]);
    % plot bayes factor time series
    plot(tif_chip, pchip(tif, log10(bf), tif_chip),  'color', 'k', 'linewidth', 1); 
    xlim(xl)
    xticks = -0.6:0.1:0;  set(gca, 'XTick', xticks);
    ylim(yl)
    yticks = -2:0.5:2;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    set(gca,'XColor','k','YColor','k')
    %set(gca, 'linewidth', 2)
    
end % condition loop

%% Row 4 (ynt sensor)
clear fslide_GA tif sem_f1 sem_f2 stat

xl = [-0.6 -0.1];

ids = 13:16;
for icond = 1:3
    axes(ha(ids(icond)));
    
    % Save data to remake figures later (for paper)
    load(fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\fslide\ynt\criterion\no_regress', ...
        sprintf('Bayes_factors_PO4_O2_PO8_%s.mat', condvect{icond})), 'bf');
    
    tif = linspace(-0.6016, -0.1211, length(bf));
    tif_chip = tif(1):0.01:tif(end);
    % add jeffrey's substantial evidence lines
    plot(tif, repmat(-0.5, size(tif)), 'color', [0.75 0.25 0.75]); hold on
    plot(tif, repmat(+0.5, size(tif)), 'color', [0.75 0.25 0.75]);
    % plot bayes factor time series
    plot(tif_chip, pchip(tif, log10(bf), tif_chip),  'color', 'k', 'linewidth', 1); 
    xlim(xl)
    xticks = -0.6:0.1:0;  set(gca, 'XTick', xticks);
    ylim(yl)
    yticks = -2:0.5:2;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    
end % condition loop


export_fig(fh2, 'D:\dfi_experiment_figures\Paper_figures\iAF\iAF_within\sd_params\fslide_criterion_no_regress_bf', '-eps', '-tiff', '-m2.5');
plot2svg('D:\dfi_experiment_figures\Paper_figures\iAF\iAF_within\sd_params\fslide_criterion_no_regress_bf_svg',fh2)












% //eof

























