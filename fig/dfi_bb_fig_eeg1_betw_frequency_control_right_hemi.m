% Figure script
% Between subject analyses
%
% Parent script(s): 
%
%   [alpha peak frequency estimates]
%       dfi_iAF_fits_eyes_open_ynt_right_hemi.m
%       dfi_iAF_fits_eyes_open_right_hemi.m
%       dfi_iAF_fits_eyes_closed_right_hemi.m
%
%   [psychometric function parameter estimates]
%       dfi_beta_binom_combine_joined_and_sep_fits.m
%   
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
% Create correlation plots between individual alpha peak 
% frequency estimates (from eyes-closed, eyes-open sensor and
% lcmv source data) and temporal window estimates (from 
% psychometric function inflection points or staircase procedure 
% in case of yes-no threshold task).
%
% Figure1: Eyes-closed
% Figure2: Eyes-open sensor level
% 
% ---
% Steffen Buergers, sbuergers@gmail.com,
% Last modified Feb. 2021


%% Setup

clear all
close all
clc

% experiment script folder
try
    addpath(genpath('dfi'))
catch
    warning('Cannot find dfi folder')
end

% run startup function
dfi_startup

% folders
data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
src_dir = fullfile(data_dir, 'source_analysis');
fig_dir  = fullfile('dfi_experiment_figures');
save_dir = fullfile(fig_dir, 'Paper_figures', 'iAF', 'iAF_betw_betabinom', ...
    'control_analysis_right_hemi'); mkdir(save_dir);
beh_figdir = fullfile('dfi_experiment_figures', 'PFs', 'beta_binom_weibull');


% add fieldtrip folder to search path
try
    addpath(fullfile('fieldtrip-20160816'))
catch
    warning('Cannot find fieldtrip folder')
end

% useful variables
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N = length(subjvect);


%% Get beta binomial data

% 2ifc
folder = '2ifc';
load(fullfile(beh_figdir, folder, 'dataPF_joined_and_sep_fits_combined.mat'), ...
    'threshold_matrix', 'slope_matrix', 'guess_matrix', 'lapse_matrix', 'eta_matrix')
threshold_2ifc = threshold_matrix;
slope_2ifc = slope_matrix;

% yesno
folder = 'yn_pooled';
load(fullfile(beh_figdir, folder, 'dataPF_joined_and_sep_fits_combined.mat'), ...
    'threshold_matrix', 'slope_matrix', 'guess_matrix', 'lapse_matrix', 'eta_matrix')
threshold_ynpool = threshold_matrix;
slope_ynpool = slope_matrix;


% Take a look at the data

disp('threshold ynpool')
min(threshold_ynpool)
max(threshold_ynpool)

disp('threshold 2ifc')
min(threshold_2ifc)
max(threshold_2ifc)


%% Load EEG data

% Corcoran eyes-closed fits
load(fullfile('dfi_experiment_figures', 'iAF_fits_corcoran_zeropadded_right_hemi', ...
    'eyes_closed_pkinfo_yn_plus_ynt.mat'))
peak_mat_ec = muPaf; clear muPaf


% Corcoran eyes-open fits (pool over yes-no and ynt)

% Load iAF peak fits (using the toolbox written by Corcoran, 2017)
load(fullfile('dfi_experiment_figures', 'iAF_fits_corcoran_zeropadded_right_hemi', ...
    'eyes_open_pkinfo.mat'));
yn_pSpec = pSpec; clear pSpec
load(fullfile('dfi_experiment_figures', 'iAF_fits_corcoran_zeropadded_right_hemi', 'yn_threshold', ...
    'eyes_open_pkinfo.mat'));
ynt_pSpec = pSpec; clear pSpec

for isubj = 1:20
    % how many sessions do we have for yn for this subject? The
    % last session (that is not nan) is the one where we estimated
    % iAF over all trials from all sessions
    sums_temp = [yn_pSpec(isubj,:).sums];
    nsess1 = length(~isnan([sums_temp.paf])) - 1;
    pSpec(isubj,1:nsess1) = yn_pSpec(isubj,1:nsess1);
    % do the same for ynt and append to pSpec
    sums_temp = [ynt_pSpec(isubj,:).sums];
    nsess2 = length(~isnan([sums_temp.paf])) - 1;
    pSpec(isubj,nsess1+1:nsess1+nsess2) = ynt_pSpec(isubj,1:nsess2); 
end

% Take a weighted average over all sessions
% weighted average of mean IAF estimates across (j-th) recordings
[muPaf, muCog] = deal(nan(20,1));
for isubj = 1:20
    Nchans = [];
    for isess = 1:size(pSpec(isubj,:),2)
        if isstruct(pSpec(isubj,isess).chans)
            Nchans = [Nchans, size(pSpec(isubj,isess).chans,2)];
        end
    end
    [muPaf(isubj, :), muCog(isubj, :)] = meanIAF_sb([pSpec(isubj,:).sums], Nchans); 
end

peak_mat_eo = muPaf; clear muPaf


% LCMV source peak estimates
% Load iAF peak fits (using the toolbox written by Corcoran, 2017)
filename = 'eyes_open_pkinfo.mat';
load(fullfile(src_dir, 'yesno', filename));
yn_pSpec = pSpec; clear pSpec
load(fullfile(src_dir, 'yn_threshold', filename));
ynt_pSpec = pSpec; clear pSpec

for isubj = 1:20
    % how many sessions do we have for yn for this subject? The
    % last session (that is not nan) is the one where we estimated
    % iAF over all trials from all sessions
    sums_temp = [yn_pSpec(isubj,:).sums];
    nsess1 = length(~isnan([sums_temp.paf])) - 1;
    pSpec(isubj,1:nsess1) = yn_pSpec(isubj,1:nsess1);
    % do the same for ynt and append to pSpec
    sums_temp = [ynt_pSpec(isubj,:).sums];
    nsess2 = length(~isnan([sums_temp.paf])) - 1;
    pSpec(isubj,nsess1+1:nsess1+nsess2) = ynt_pSpec(isubj,1:nsess2); 
end

% Take a weighted average over all sessions
% weighted average of mean IAF estimates across (j-th) recordings
[muPaf, muCog] = deal(nan(20,1));
for isubj = 1:20
    Nchans = [];
    for isess = 1:size(pSpec(isubj,:),2)
        if isstruct(pSpec(isubj,isess).chans)
            Nchans = [Nchans, size(pSpec(isubj,isess).chans,2)];
        end
    end
    [muPaf(isubj, :), muCog(isubj, :)] = meanIAF_sb([pSpec(isubj,:).sums], Nchans); 
end
peak_mat_occ_source = muPaf; clear muPaf


% Some settings
col_vect = [[0 0.6 0]; [0 0 1]; [1 0 0]; [0 0 0]];
col_vect_ci = [[0.75 1 0.75]; [0.75 0.75 1]; [1 0.75 0.75]; [0.75 0.75 0.75]];
col_lines = [[0 0.8 0]; [0 0 1]; [1 0 0]; [0 0 0]];
opacity = 0.5; % transparency (alpha)
lw = 0.6;
ls = '-';


%% Figure 0: Threshold: Separate scatter plots for conditions, only sensory eyes-open

% Columns: iAF estimates
% Rows: 2IFC, yes-no pooled, staircase SOA

fh0 = figure('color', [1 1 1], 'Position', [0, 0, 427, 350]);
ha = tight_subplot(3, 4, [0.01 0.03], [0.02], [0.02]);

plot_ci = false;

% 1.) 2IFC threshold vs iAF

[spearRho, pval, nobs, r, b1, b0] = deal(zeros(1, 3));
for icond = 1:3

    % [3,1] Threshold 2ifc vs eyes-closed Corcoran iAF
    ip = threshold_2ifc;
    yl = [0, 0.245]; yrange = yl(2) - yl(1);
    
    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(ip) | repmat(~any(~isnan(peak_mat_eo),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);
        
    ip(~good_fit) = nan;

    % Correlate iAF with perceptual windows for all tasks (most importantly dfi)
    percwin = ones(size(peak_mat_eo))./peak_mat_eo;

    axes(ha(icond));
    percwin_mat = repmat(percwin, [1,3]);
    percwin_mat(~good_fit) = nan;
    ip(~good_fit) = nan;
    xl = [min(min(1./percwin_mat))*0.95, max(max(1./percwin_mat))*1.05]; xrange = xl(2) - xl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(1./percwin_mat(:,icond), ip(:,icond), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(1./percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));
        line([min(min(1./percwin_mat))  max(max(1./percwin_mat))], [b0(1,ic)+b1(1,ic)*min(min(1./percwin_mat)) b0(1,ic)+b1(1,ic)*max(max(1./percwin_mat))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = 6:1:14; set(gca, 'XTick', xticks)
        yticks = -0.1:0.1:0.2;  set(gca, 'YTick', yticks);
        if plot_ci
            % Obtain prediction intervals, and plot them
            [P,S] = polyfit(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic),1);
            xvals = [min(min(1./percwin_mat)):0.001:max(max(1./percwin_mat))];
            [Y,DELTA] = polyconf([b1(1,ic), b0(1,ic)],xvals,S);
            plot(xvals, Y+DELTA./1000, '--', 'color', col_vect_ci(ic,:)); hold on
            plot(xvals, Y-DELTA./1000, '--', 'color', col_vect_ci(ic,:))
        end
    end
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])

end

cndttl = {'0S', '1S', '2S'};
disp('------------------------------------')
disp('eyes-open --- 2IFC threshold vs iAF:')
disp('------------------------------------')
for ic = 1:3
    fprintf('%s   -->   N = %i,     r = %f,      BF = %f\n', ...
        cndttl{ic}, nobs(1,ic), r(1,ic), corrbf(r(1,ic), nobs(1,ic)))
end
fprintf('\n\n\n')

% add Bayes factors
axes(ha(4))
ic=1; line([0 0],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0.6 0], 'linewidth', 6); hold on
ic=2; line([0.5 0.5],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0 1], 'linewidth', 6)
ic=3; line([1 1],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[1 0 0], 'linewidth', 6)
line([-0.5 1.5], [0 0], 'color', 'k')
line([-0.5 1.5], [0.5 0.5], 'color', [0.75 0.25 0.75])
line([-0.5 1.5], [-0.5 -0.5], 'color', [0.75 0.25 0.75])
xlim([-0.5 2.5])
ylim([-1.2 1.2])


% 2.) yes-no pooled threshold vs iAF

[spearRho, pval, nobs, r, b1, b0] = deal(zeros(1, 3));
for icond = 1:3

    % [3,1] Threshold 2ifc vs eyes-closed Corcoran iAF
    ip = threshold_ynpool;
    yl = [0, 0.245]; yrange = yl(2) - yl(1);
    
    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(ip) | repmat(~any(~isnan(peak_mat_eo),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);
        
    ip(~good_fit) = nan;

    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(ip) | repmat(~any(~isnan(peak_mat_eo),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);

    % Correlate iAF with perceptual windows for all tasks (most importantly dfi)
    percwin = ones(size(peak_mat_eo))./peak_mat_eo;

    axes(ha(icond+4));
    percwin_mat = repmat(percwin, [1,3]);
    percwin_mat(~good_fit) = nan;
    ip(~good_fit) = nan;
    xl = [min(min(1./percwin_mat))*0.95, max(max(1./percwin_mat))*1.05]; xrange = xl(2) - xl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(1./percwin_mat(:,icond), ip(:,icond), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(1./percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));
        line([min(min(1./percwin_mat))  max(max(1./percwin_mat))], [b0(1,ic)+b1(1,ic)*min(min(1./percwin_mat)) b0(1,ic)+b1(1,ic)*max(max(1./percwin_mat))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = 6:1:14; set(gca, 'XTick', xticks)
        yticks = -0.1:0.1:0.2;  set(gca, 'YTick', yticks);
        if plot_ci
            % Obtain prediction intervals, and plot them
            [P,S] = polyfit(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic),1);
            xvals = [min(min(1./percwin_mat)):0.001:max(max(1./percwin_mat))];
            [Y,DELTA] = polyconf([b1(1,ic), b0(1,ic)],xvals,S);
            plot(xvals, Y+DELTA./1000, '--', 'color', col_vect_ci(ic,:)); hold on
            plot(xvals, Y-DELTA./1000, '--', 'color', col_vect_ci(ic,:))
        end
    end
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])

end

cndttl = {'0S', '1S', '2S'};
disp('---------------------------------------------')
disp('eyes-open --- yes-no pooled threshold vs iAF:')
disp('---------------------------------------------')
for ic = 1:3
    fprintf('%s   -->   N = %i,     r = %f,      BF = %f\n', ...
        cndttl{ic}, nobs(1,ic), r(1,ic), corrbf(r(1,ic), nobs(1,ic)))
end
fprintf('\n\n\n')


% add Bayes factors
axes(ha(4+4))
ic=1; line([0 0],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0.6 0], 'linewidth', 6); hold on
ic=2; line([0.5 0.5],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0 1], 'linewidth', 6)
ic=3; line([1 1],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[1 0 0], 'linewidth', 6)
line([-0.5 1.5], [0 0], 'color', 'k')
line([-0.5 1.5], [0.5 0.5], 'color', [0.75 0.25 0.75])
line([-0.5 1.5], [-0.5 -0.5], 'color', [0.75 0.25 0.75])
xlim([-0.5 2.5])
ylim([-1.2 1.2])


% 3.) yn-threshold SOA versus iAF

load(fullfile('dfi_experiment_data', 'data', 'experiment', 'd701to727_ynt.mat'))
clear d7*

subjects = unique(dall.partid);
trialtypes = [3, 6, 8, 9];
soamat = nan(20,4);
for isubj = 1:N
    for itrl = 1:length(trialtypes)
        soamat(isubj,itrl) = unique(dall.soa(dall.partid == subjects(isubj) & dall.trlid == trialtypes(itrl)) );
    end
end

[spearRho, pval, nobs, r, b1, b0] = deal(zeros(1, 3));
for icond = 1:3

    % [1] Threshold yes-no pool vs eyes-open Corcoran iAF
    % dfi and dfi_control have the same soa anyway
    ip = soamat(:,1:3);
    %ip = threshold_ynpool;
    yl = [0, 0.245]; yrange = yl(2) - yl(1);

    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(ip) | repmat(~any(~isnan(peak_mat_eo),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);

    % Correlate iAF with perceptual windows for all tasks (most importantly dfi)
    percwin = ones(size(peak_mat_eo))./peak_mat_eo;

    axes(ha(icond+8));
    percwin_mat = repmat(percwin, [1,3]);
    percwin_mat(~good_fit) = nan;
    ip(~good_fit) = nan;
    xl = [min(min(1./percwin_mat))*0.95, max(max(1./percwin_mat))*1.05]; xrange = xl(2) - xl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(1./percwin_mat(:,icond), ip(:,icond), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(1./percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));
        line([min(min(1./percwin_mat))  max(max(1./percwin_mat))], [b0(1,ic)+b1(1,ic)*min(min(1./percwin_mat)) b0(1,ic)+b1(1,ic)*max(max(1./percwin_mat))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = 6:1:14; set(gca, 'XTick', xticks)
        yticks = -0.1:0.1:0.2;  set(gca, 'YTick', yticks);
        if plot_ci
            % Obtain prediction intervals, and plot them
            [P,S] = polyfit(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic),1);
            xvals = [min(min(1./percwin_mat)):0.001:max(max(1./percwin_mat))];
            [Y,DELTA] = polyconf([b1(1,ic), b0(1,ic)],xvals,S);
            plot(xvals, Y+DELTA./1000, '--', 'color', col_vect_ci(ic,:)); hold on
            plot(xvals, Y-DELTA./1000, '--', 'color', col_vect_ci(ic,:))
        end
    end
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])

end

cndttl = {'0S', '1S', '2S'};
disp('------------------------------------------')
disp('eyes-open --- yn-threshold SOA versus iAF:')
disp('------------------------------------------')
for ic = 1:3
    fprintf('%s   -->   N = %i,     r = %f,      BF = %f\n', ...
        cndttl{ic}, nobs(1,ic), r(1,ic), corrbf(r(1,ic), nobs(1,ic)))
end
fprintf('\n\n\n')

% add Bayes factors
axes(ha(4+8))
ic=1; line([0 0],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0.6 0], 'linewidth', 6); hold on
ic=2; line([0.5 0.5],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0 1], 'linewidth', 6)
ic=3; line([1 1],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[1 0 0], 'linewidth', 6)
line([-0.5 1.5], [0 0], 'color', 'k')
line([-0.5 1.5], [0.5 0.5], 'color', [0.75 0.25 0.75])
line([-0.5 1.5], [-0.5 -0.5], 'color', [0.75 0.25 0.75])
xlim([-0.5 2.5])
ylim([-1.2 1.2])


box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.02])

saveas(fh0,fullfile(save_dir, 'eyes_open_sensor_threshold_scatter_svg_freq.svg'))

close all


%% Figure 01: Threshold: Separate scatter plots for conditions, only sensory eyes-closed

% Columns: iAF estimates
% Rows: 2IFC, yes-no pooled, staircase SOA

fh01 = figure('color', [1 1 1], 'Position', [0, 0, 427, 350]);
ha = tight_subplot(3, 4, [0.01 0.03], [0.02], [0.02]);

plot_ci = false;

% 1.) 2IFC threshold vs iAF

[spearRho, pval, nobs, r, b1, b0] = deal(zeros(1, 3));
for icond = 1:3

    % [3,1] Threshold 2ifc vs eyes-closed Corcoran iAF
    ip = threshold_2ifc;
    yl = [0, 0.245]; yrange = yl(2) - yl(1);
    
    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(ip) | repmat(~any(~isnan(peak_mat_ec),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);
        
    ip(~good_fit) = nan;

    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(ip) | repmat(~any(~isnan(peak_mat_ec),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);

    % Correlate iAF with perceptual windows for all tasks (most importantly dfi)
    percwin = ones(size(peak_mat_ec))./peak_mat_ec;

    axes(ha(icond));
    percwin_mat = repmat(percwin, [1,3]);
    percwin_mat(~good_fit) = nan;
    ip(~good_fit) = nan;
    xl = [min(min(1./percwin_mat))*0.95, max(max(1./percwin_mat))*1.05]; xrange = xl(2) - xl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(1./percwin_mat(:,icond), ip(:,icond), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(1./percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));
        line([min(min(1./percwin_mat))  max(max(1./percwin_mat))], [b0(1,ic)+b1(1,ic)*min(min(1./percwin_mat)) b0(1,ic)+b1(1,ic)*max(max(1./percwin_mat))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = 6:1:14; set(gca, 'XTick', xticks)
        yticks = -0.1:0.1:0.2;  set(gca, 'YTick', yticks);
        if plot_ci
            % Obtain prediction intervals, and plot them
            [P,S] = polyfit(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic),1);
            xvals = [min(min(1./percwin_mat)):0.001:max(max(1./percwin_mat))];
            [Y,DELTA] = polyconf([b1(1,ic), b0(1,ic)],xvals,S);
            plot(xvals, Y+DELTA./1000, '--', 'color', col_vect_ci(ic,:)); hold on
            plot(xvals, Y-DELTA./1000, '--', 'color', col_vect_ci(ic,:))
        end
    end
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])

end

cndttl = {'0S', '1S', '2S'};
disp('------------------------------------------')
disp('eyes-closed --- 2IFC threshold vs iAF:')
disp('------------------------------------------')
for ic = 1:3
    fprintf('%s   -->   N = %i,     r = %f,      BF = %f\n', ...
        cndttl{ic}, nobs(1,ic), r(1,ic), corrbf(r(1,ic), nobs(1,ic)))
end
fprintf('\n\n\n')

% add Bayes factors
axes(ha(4))
ic=1; line([0 0],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0.6 0], 'linewidth', 6); hold on
ic=2; line([0.5 0.5],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0 1], 'linewidth', 6)
ic=3; line([1 1],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[1 0 0], 'linewidth', 6)
line([-0.5 1.5], [0 0], 'color', 'k')
line([-0.5 1.5], [0.5 0.5], 'color', [0.75 0.25 0.75])
line([-0.5 1.5], [-0.5 -0.5], 'color', [0.75 0.25 0.75])
xlim([-0.5 2.5])
ylim([-1.2 1.2])

box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.02])


% 2.) yes-no pooled threshold vs iAF

[spearRho, pval, nobs, r, b1, b0] = deal(zeros(1, 3));
for icond = 1:3

    % [3,1] Threshold 2ifc vs eyes-closed Corcoran iAF
    ip = threshold_ynpool;
    yl = [0, 0.245]; yrange = yl(2) - yl(1);
    
    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(ip) | repmat(~any(~isnan(peak_mat_ec),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);
        
    ip(~good_fit) = nan;

    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(ip) | repmat(~any(~isnan(peak_mat_ec),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);

    % Correlate iAF with perceptual windows for all tasks (most importantly dfi)
    percwin = ones(size(peak_mat_ec))./peak_mat_ec;

    axes(ha(icond+4));
    percwin_mat = repmat(percwin, [1,3]);
    percwin_mat(~good_fit) = nan;
    ip(~good_fit) = nan;
    xl = [min(min(1./percwin_mat))*0.95, max(max(1./percwin_mat))*1.05]; xrange = xl(2) - xl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(1./percwin_mat(:,icond), ip(:,icond), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(1./percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));
        line([min(min(1./percwin_mat))  max(max(1./percwin_mat))], [b0(1,ic)+b1(1,ic)*min(min(1./percwin_mat)) b0(1,ic)+b1(1,ic)*max(max(1./percwin_mat))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = 6:1:14; set(gca, 'XTick', xticks)
        yticks = -0.1:0.1:0.2;  set(gca, 'YTick', yticks);
        if plot_ci
            % Obtain prediction intervals, and plot them
            [P,S] = polyfit(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic),1);
            xvals = [min(min(1./percwin_mat)):0.001:max(max(1./percwin_mat))];
            [Y,DELTA] = polyconf([b1(1,ic), b0(1,ic)],xvals,S);
            plot(xvals, Y+DELTA./1000, '--', 'color', col_vect_ci(ic,:)); hold on
            plot(xvals, Y-DELTA./1000, '--', 'color', col_vect_ci(ic,:))
        end
    end
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])

end

cndttl = {'0S', '1S', '2S'};
disp('-----------------------------------------------')
disp('eyes-closed --- yes-no pooled threshold vs iAF:')
disp('-----------------------------------------------')
for ic = 1:3
    fprintf('%s   -->   N = %i,     r = %f,      BF = %f\n', ...
        cndttl{ic}, nobs(1,ic), r(1,ic), corrbf(r(1,ic), nobs(1,ic)))
end
fprintf('\n\n\n')

% add Bayes factors
axes(ha(4+4))
ic=1; line([0 0],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0.6 0], 'linewidth', 6); hold on
ic=2; line([0.5 0.5],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0 1], 'linewidth', 6)
ic=3; line([1 1],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[1 0 0], 'linewidth', 6)
line([-0.5 1.5], [0 0], 'color', 'k')
line([-0.5 1.5], [0.5 0.5], 'color', [0.75 0.25 0.75])
line([-0.5 1.5], [-0.5 -0.5], 'color', [0.75 0.25 0.75])
xlim([-0.5 2.5])
ylim([-1.2 1.2])

box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.02])


% 3.) yn-threshold SOA versus iAF

load(fullfile('dfi_experiment_data', 'data', 'experiment', 'd701to727_ynt.mat'))
clear d7*

subjects = unique(dall.partid);
trialtypes = [3, 6, 8, 9];
soamat = nan(20,4);
for isubj = 1:N
    for itrl = 1:length(trialtypes)
        soamat(isubj,itrl) = unique(dall.soa(dall.partid == subjects(isubj) & dall.trlid == trialtypes(itrl)) );
    end
end

[spearRho, pval, nobs, r, b1, b0] = deal(zeros(1, 3));
for icond = 1:3

    % [1] Threshold yes-no pool vs eyes-open Corcoran iAF
    % dfi and dfi_control have the same soa anyway
    ip = soamat(:,1:3);
    %ip = threshold_ynpool;
    yl = [0, 0.245]; yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 

    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(ip) | repmat(~any(~isnan(peak_mat_ec),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);

    % Correlate iAF with perceptual windows for all tasks (most importantly dfi)
    percwin = ones(size(peak_mat_ec))./peak_mat_ec;

    axes(ha(icond+8));
    percwin_mat = repmat(percwin, [1,3]);
    percwin_mat(~good_fit) = nan;
    ip(~good_fit) = nan;
    xl = [min(min(1./percwin_mat))*0.95, max(max(1./percwin_mat))*1.05]; xrange = xl(2) - xl(1);
    transparentScatter(1./percwin_mat(:,icond), ip(:,icond), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(1./percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));
        line([min(min(1./percwin_mat))  max(max(1./percwin_mat))], [b0(1,ic)+b1(1,ic)*min(min(1./percwin_mat)) b0(1,ic)+b1(1,ic)*max(max(1./percwin_mat))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = 6:1:14; set(gca, 'XTick', xticks)
        yticks = -0.1:0.1:0.2;  set(gca, 'YTick', yticks);
        if plot_ci
            % Obtain prediction intervals, and plot them
            [P,S] = polyfit(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic),1);
            xvals = [min(min(1./percwin_mat)):0.001:max(max(1./percwin_mat))];
            [Y,DELTA] = polyconf([b1(1,ic), b0(1,ic)],xvals,S);
            plot(xvals, Y+DELTA./1000, '--', 'color', col_vect_ci(ic,:)); hold on
            plot(xvals, Y-DELTA./1000, '--', 'color', col_vect_ci(ic,:))
        end
    end
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])

end

cndttl = {'0S', '1S', '2S'};
disp('--------------------------------------------')
disp('eyes-closed --- yn-threshold SOA versus iAF:')
disp('--------------------------------------------')
for ic = 1:3
    fprintf('%s   -->   N = %i,     r = %f,      BF = %f\n', ...
        cndttl{ic}, nobs(1,ic), r(1,ic), corrbf(r(1,ic), nobs(1,ic)))
end
fprintf('\n\n\n')

% add Bayes factors
axes(ha(4+8))
ic=1; line([0 0],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0.6 0], 'linewidth', 6); hold on
ic=2; line([0.5 0.5],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[0 0 1], 'linewidth', 6)
ic=3; line([1 1],[0 log10(corrbf(r(1,ic), nobs(1,ic)))],'color',[1 0 0], 'linewidth', 6)
line([-0.5 1.5], [0 0], 'color', 'k')
line([-0.5 1.5], [0.5 0.5], 'color', [0.75 0.25 0.75])
line([-0.5 1.5], [-0.5 -0.5], 'color', [0.75 0.25 0.75])
xlim([-0.5 2.5])
ylim([-1.2 1.2])

box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.02])

saveas(fh01,fullfile(save_dir, 'eyes_closed_sensor_threshold_scatter_svg_freq.svg'))

close all


% eof

