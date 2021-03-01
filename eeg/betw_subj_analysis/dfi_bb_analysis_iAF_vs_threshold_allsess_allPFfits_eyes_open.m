% Assess relationship between individual alpha frequency from 
% task recordings and psychometric function inflection points.
%
% Parent script(s): 
%   dfi_iAF_fits_eyes_open.m
%   dfi_iAF_fits_eyes_open_ynt.m
%   dfi_compare_betabinom_and_binom_scripts.m
%
% Children script(s): 
%   ?
%
% Sibling script(s):
%   dfi_bb_fig_eeg1_betw_frequency.m
%
% DETAILS
%
%   ?
%
% NOTES
% Be aware that here we compute correlations between wavelengths and 
% inflection points of psychometric functions. In
% dfi_bb_fig_eeg1_betw_frequency.m we use frequency instead. This yields
% slightly different statistics. For the manuscript we use frequency, not
% wavelength!
% 
% ===========================================================================
%
%     dfi (double flash illusion) codebase accompanying the manuscript ...
%     Copyright (C) 2021  Steffen Buergers
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% ---
% Steffen Buergers, sbuergers@gmail.com,
% Last modified Feb. 2021


%% *** SETUP ***
clear all

% disable warnings for speed sakD:
% warning('off','all')

% Name of directory to save things to
an_fold = 'iAF_percwin_between_zeropadded_allPFfits_beta_binom';

% save figures?
save_figures = true;

% experiment script folder
try
    addpath(genpath('dfi'))
catch
    warning('Cannot find dfi folder')
end

% run startup function
dfi_startup

% experiment data folder
data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
fig_dir  = fullfile('dfi_experiment_figures');

% beta binom data dir
figdir = fullfile('dfi_experiment_figures', 'PFs', 'beta_binom_weibull');

% add fieldtrip folder to search path
try
    addpath(fullfile('toolboxes', 'fieldtrip-20160816'))
catch
    warning('Cannot find fieldtrip folder')
end

% useful variables
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
filename = 'peaks_corcoran'; % peaks_search2peaks
N = length(subjvect);
          
% Load iAF peak fits (using the toolbox written by Corcoran, 2017)
%load(fullfile(fig_dir, 'iAF_fits_corcoran', 'eyes_open_pkinfo.mat'));

% Load iAF peak fits (using the toolbox written by Corcoran, 2017)
load(fullfile('dfi_experiment_figures', 'iAF_fits_corcoran_zeropadded', 'eyes_open_pkinfo.mat'));
yn_pSpec = pSpec; clear pSpec
load(fullfile('dfi_experiment_figures', 'iAF_fits_corcoran_zeropadded', 'yn_threshold', 'eyes_open_pkinfo.mat'));
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

% We might look at single sessions, or all of them combined (for now I only
% look at all of them combined, because I don't have psychometric function
% fits for single sessions...)
nfiles = 1;
for ifold = 1:nfiles
    
    % Prepare folders and peak matrices
    if ifold == nfiles
        fold_data  = 'all_yn_sessions';
        peak_mat = muPaf;
    else
        % ...
    end
    
    % make folders
    if ~exist(fullfile(fig_dir, an_fold, fold_data), 'dir')
        mkdir(fullfile(fig_dir, an_fold, fold_data));
    end

    
    %% *** YN-pooled (logistic) THRESHOLD ***
    
    % Get psychometric function data
    % Here I already excluded bad fits. Plus separate function fits were
    % only substituted where the overall/joined fits failed!
    % Import beta binomial parameter estimates
    folder = 'yn_pooled';
    load(fullfile(figdir, folder, 'data_betabinom_in_pffit_format_yesno.mat'));
    pffit_all = pffit_yn;
    
    ip = nan(length(pffit_all), length(pffit_all{1}));
    for isubj = 1:N
        for icond = 1:3
            if ~isempty(pffit_all{isubj}{icond})
                ip(isubj,icond) = pffit_all{isubj}{icond}.par(1);
            end
        end
    end
      
    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(ip) | repmat(~any(~isnan(peak_mat),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);
    
    % Plot correlations between perceptual window and iAF per channel and
    % condition
    percwin = ones(size(peak_mat))./peak_mat;
    col_vect = [[0.4 1 0.4]; [0.4 0.4 1]; [1 0.4 0.4]];
    [spearRho, pval, nobs, r, b1, b0] = deal(zeros(64,3));
    fh = figure('color', 'w', 'position', [50 50 800, 790]);
    for ic = 1:3
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));

        subplot(2,2,ic)
        plot(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), '.', 'color', col_vect(ic,:), 'markersize', 20); hold on
        line([1/12  1/7], [b0(1,ic)+b1(1,ic)*1/12 b0(1,ic)+b1(1,ic)*1/7], 'color', col_vect(ic,:)./2);
        xlim([1/12, 1/7]); ylim([-0.125 0.25]);
        xticks = 0.07:0.01:0.14; set(gca, 'XTick', xticks)
        yticks = -0.1:0.05:0.2; set(gca, 'YTick', yticks);
        addtext(sprintf('\nr^2 = %.2f\niAF = %.2f + %.2f*pWin\np-value = %.2f\nN = %.2f', spearRho(1,ic), b0(1,ic), b1(1,ic), pval(1,ic), nobs(1,ic)));
        if ic == 1
            title('1 Flash & 2 Flashes');
        elseif ic == 2
            title('1 Flash, 1 Sound & 2 Flashes, 1 Sound');
        elseif ic == 3
            title('1 Flash, 2 Sounds & 2 Flashes, 2 Sounds');
        end
        if ic > 2, xlabel('Occipital alpha cycle duration (s)'); end
        if ~isEven(ic), ylabel('Inflection point (s)'); end
        plotspecs
    end

    if ~exist(fullfile(fig_dir, an_fold, fold_data), 'dir');
        mkdir(fullfile(fig_dir, an_fold, fold_data));
    end
    saveas(fh, fullfile(fig_dir, an_fold, fold_data, sprintf('yn_pooled_iAF_perWin_corr_eyes_open.emf')))
    close all
    save(fullfile(fig_dir, an_fold, fold_data, 'iAF_perWin_data_yn_pooled_eyes_open'), 'spearRho', 'pval', 'r', 'b1', 'b0', 'nobs')
    
    
    
    %% *** YN-pooled (logistic) SLOPE ***

    % Get psychometric function data
    % Here I already excluded bad fits. Plus separate function fits were
    % only substituted where the overall/joined fits failed!
    
    slope = nan(length(pffit_all), length(pffit_all{1}));
    for isubj = 1:N
        for icond = 1:3
            if ~isempty(pffit_all{isubj}{icond})
                slope(isubj,icond) = pffit_all{isubj}{icond}.par(2);
            end
        end
    end
    
    % Plot correlations between perceptual window and iAF per channel and
    % condition
    percwin = ones(size(peak_mat))./peak_mat;
    col_vect = [[0.4 1 0.4]; [0.4 0.4 1]; [1 0.4 0.4]];
    [spearRho, pval, nobs, r, b1, b0] = deal(zeros(64,3));
    fh = figure('color', 'w', 'position', [50 50 800, 790]);
    for ic = 1:3
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));

        subplot(2,2,ic)
        plot(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), '.', 'color', col_vect(ic,:), 'markersize', 20); hold on
        line([1/12  1/7], [b0(1,ic)+b1(1,ic)*1/12 b0(1,ic)+b1(1,ic)*1/7], 'color', col_vect(ic,:)./2);
        xlim([1/12, 1/7]); ylim([0 5]);
        xticks = 0.07:0.01:0.14; set(gca, 'XTick', xticks)
        yticks = -5:0.5:5; set(gca, 'YTick', yticks);
        addtext(sprintf('\nr^2 = %.2f\niAF = %.2f + %.2f*pWin\np-value = %.2f\nN = %.2f', spearRho(1,ic), b0(1,ic), b1(1,ic), pval(1,ic), nobs(1,ic)));
        if ic == 1
            title('1 Flash & 2 Flashes');
        elseif ic == 2
            title('1 Flash, 1 Sound & 2 Flashes, 1 Sound');
        elseif ic == 3
            title('1 Flash, 2 Sounds & 2 Flashes, 2 Sounds');
        end
        if ic > 2, xlabel('Occipital alpha cycle duration (s)'); end
        if ~isEven(ic), ylabel('log-slope'); end
        plotspecs
    end

    if ~exist(fullfile(fig_dir, an_fold, fold_data), 'dir');
        mkdir(fullfile(fig_dir, an_fold, fold_data));
    end
    saveas(fh, fullfile(fig_dir, an_fold, fold_data, sprintf('yn_pooled_iAF_slope_corr_eyes_open.emf')))
    close all
    save(fullfile(fig_dir, an_fold, fold_data, 'iAF_slope_data_yn_pooled_eyes_open'), 'spearRho', 'pval', 'r', 'b1', 'b0', 'nobs')
    
    
    
    
    
    
    %% *** 2IFC (logistic) THRESHOLD ***

    % Get psychometric function data
    % Here I already excluded bad fits. Plus separate function fits were
    % only substituted where the overall/joined fits failed!
    folder = '2ifc';
    load(fullfile(figdir, folder, 'data_betabinom_in_pffit_format_2IFC.mat'));
    pffit_all = pffit_2ifc;
    
    ip = nan(length(pffit_all), length(pffit_all{1}));
    for isubj = 1:N
        for icond = 1:3
            if ~isempty(pffit_all{isubj}{icond})
                ip(isubj,icond) = pffit_all{isubj}{icond}.par(1);
            end
        end
    end
      
    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(ip) | repmat(~any(~isnan(peak_mat),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);
    
    % Plot correlations between perceptual window and iAF per channel and
    % condition
    percwin = ones(size(peak_mat))./peak_mat;
    col_vect = [[0.4 1 0.4]; [0.4 0.4 1]; [1 0.4 0.4]];
    [spearRho, pval, nobs, r, b1, b0] = deal(zeros(64,3));
    fh = figure('color', 'w', 'position', [50 50 800, 790]);
    for ic = 1:3
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));

        subplot(2,2,ic)
        plot(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), '.', 'color', col_vect(ic,:), 'markersize', 20); hold on
        line([1/12  1/7], [b0(1,ic)+b1(1,ic)*1/12 b0(1,ic)+b1(1,ic)*1/7], 'color', col_vect(ic,:)./2);
        xlim([1/12, 1/7]); ylim([-0.125 0.25]);
        xticks = 0.07:0.01:0.14; set(gca, 'XTick', xticks)
        yticks = -0.1:0.05:0.2; set(gca, 'YTick', yticks);
        addtext(sprintf('\nr^2 = %.2f\niAF = %.2f + %.2f*pWin\np-value = %.2f\nN = %.2f', spearRho(1,ic), b0(1,ic), b1(1,ic), pval(1,ic), nobs(1,ic)));
        if ic == 1
            title('1 Flash & 2 Flashes');
        elseif ic == 2
            title('1 Flash, 1 Sound & 2 Flashes, 1 Sound');
        elseif ic == 3
            title('1 Flash, 2 Sounds & 2 Flashes, 2 Sounds');
        end
        if ic > 2, xlabel('Occipital alpha cycle duration (s)'); end
        if ~isEven(ic), ylabel('Inflection point (s)'); end
        plotspecs
    end

    if ~exist(fullfile(fig_dir, an_fold, fold_data), 'dir');
        mkdir(fullfile(fig_dir, an_fold, fold_data));
    end
    saveas(fh, fullfile(fig_dir, an_fold, fold_data, sprintf('twoifc_iAF_perWin_corr_eyes_open.emf')))
    close all
    save(fullfile(fig_dir, an_fold, fold_data, 'iAF_perWin_data_twoifc_eyes_open'), 'spearRho', 'pval', 'r', 'b1', 'b0', 'nobs')
    

    
    
    %% *** 2IFC (logistic) SLOPE ***

    % Get psychometric function data
    % Here I already excluded bad fits. Plus separate function fits were
    % only substituted where the overall/joined fits failed!
    slope = nan(length(pffit_all), length(pffit_all{1}));
    for isubj = 1:N
        for icond = 1:3
            if ~isempty(pffit_all{isubj}{icond})
                slope(isubj,icond) = pffit_all{isubj}{icond}.par(2);
            end
        end
    end
    
    % Plot correlations between perceptual window and iAF per channel and
    % condition
    percwin = ones(size(peak_mat))./peak_mat;
    col_vect = [[0.4 1 0.4]; [0.4 0.4 1]; [1 0.4 0.4]];
    [spearRho, pval, nobs, r, b1, b0] = deal(zeros(64,3));
    fh = figure('color', 'w', 'position', [50 50 800, 790]);
    for ic = 1:3
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));

        subplot(2,2,ic)
        plot(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), '.', 'color', col_vect(ic,:), 'markersize', 20); hold on
        line([1/12  1/7], [b0(1,ic)+b1(1,ic)*1/12 b0(1,ic)+b1(1,ic)*1/7], 'color', col_vect(ic,:)./2);
        xlim([1/12, 1/7]); ylim([0 5]);
        xticks = 0.07:0.01:0.14; set(gca, 'XTick', xticks)
        yticks = -5:0.5:5; set(gca, 'YTick', yticks);
        addtext(sprintf('\nr^2 = %.2f\niAF = %.2f + %.2f*pWin\np-value = %.2f\nN = %.2f', spearRho(1,ic), b0(1,ic), b1(1,ic), pval(1,ic), nobs(1,ic)));
        if ic == 1
            title('1 Flash & 2 Flashes');
        elseif ic == 2
            title('1 Flash, 1 Sound & 2 Flashes, 1 Sound');
        elseif ic == 3
            title('1 Flash, 2 Sounds & 2 Flashes, 2 Sounds');
        end
        if ic > 2, xlabel('Occipital alpha cycle duration (s)'); end
        if ~isEven(ic), ylabel('log-slope'); end
        plotspecs
    end

    if ~exist(fullfile(fig_dir, an_fold, fold_data), 'dir');
        mkdir(fullfile(fig_dir, an_fold, fold_data));
    end
    saveas(fh, fullfile(fig_dir, an_fold, fold_data, sprintf('twoifc_iAF_slope_corr_eyes_open.emf')))
    close all
    save(fullfile(fig_dir, an_fold, fold_data, 'iAF_slope_data_twoifc_eyes_open'), 'spearRho', 'pval', 'r', 'b1', 'b0', 'nobs')
    
    
    
    %% *** YNT data. SOA of PSE ***

    % Get psychometric function data
    load('D:\dfi_experiment_data\data\experiment\d701to727_ynt.mat')
    clear d7*
    
    subjects = unique(dall.partid);
    trialtypes = [3, 6, 8, 9];
    soamat = nan(20,4);
    for isubj = 1:N
        for itrl = 1:length(trialtypes)
            soamat(isubj,itrl) = unique(dall.soa(dall.partid == subjects(isubj) & dall.trlid == trialtypes(itrl)) );
        end
    end
    
    figure;
    plot(1:20, soamat(:,1), 'g'); hold on
    plot(1:20, soamat(:,2), 'b');
    plot(1:20, soamat(:,3), 'r');
    
    
    % For now, exclude non-sensical / clearly bad fits
    bad_fit  = isnan(peak_mat);
    good_fit = repmat(~(bad_fit), [1,4]);
    
    ip = soamat;
    
    % Correlate iAF with perceptual windows for all tasks (most importantly dfi)
    percwin = ones(size(peak_mat))./peak_mat;
    fh = figure('color', 'w', 'position', [50 50 800, 790]);
    col_vect = [[0.4 1 0.4]; [0.4 0.4 1]; [1 0.4 0.4]; [0.4 0.4 0.4]];
    [pearsR, spearRho, pval, pval_r, nobs, r, b1, b0] = deal(zeros(4,1));
    for ic = 1:4
        [pearsR(ic), pval_r(ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic));
        [spearRho(ic), pval(ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(ic),b1(ic),b0(ic)] = regression(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'one');
        nobs(ic) = sum(good_fit(:,ic));
        
        subplot(2,2,ic)
        plot(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), '.', 'color', col_vect(ic,:), 'markersize', 20); hold on
        line([1/14 1/8], [b0(ic)+b1(ic)*1/14 b0(ic)+b1(ic)*1/8], 'color', col_vect(ic,:)./2);
        xlim([1/14, 1/8]); ylim([0 0.225]);
        xticks = 0.07:0.01:0.12; set(gca, 'XTick', xticks)
        yticks = 0:0.05:0.225;         set(gca, 'YTick', yticks);
        addtext(sprintf('r^2 = %.2f\niAF = %.2f + %.2f*pWin\np-value = %.2f\nN = %.2f', pearsR(ic), b0(ic), b1(ic), pval(ic), nobs(ic)));
        ylim([0 0.25]);
        if ic == 1
            title('2 Flashes');
        elseif ic == 2
            title('2 Flashes, 1 Sound');
        elseif ic == 3
            title('1 Flash, 2 Sounds');
        else
            title('2 Flashes, 2 Sounds');
        end
        if ic > 2, xlabel('Occipital alpha cycle duration (s)'); end
        if ~isEven(ic), ylabel('Inflection point (s)'); end
        plotspecs
    end
    [pearsR, pval, nobs]
    foldername = 'ynt_soa';
    if ~exist(fullfile(fig_dir, an_fold, fold_data), 'dir');
        mkdir(fullfile(fig_dir, an_fold, fold_data));
    end
    saveas(fh, fullfile(fig_dir, an_fold, fold_data, sprintf('iAF_vs_ynt_soa.emf')))
    close all
    save(fullfile(fig_dir, an_fold, fold_data, 'iAF_vs_ynt_soa'), 'spearRho', 'pval', 'r', 'b1', 'b0', 'nobs')
    
    
    
    %% *** YN data. Response times 2IFC ***
    
    % Get response time data
    load('D:\dfi_experiment_data\data\experiment\d701to727_2ifc.mat')
    dall.trlid(dall.trlid == 2) = 3;
    dall.trlid(dall.trlid == 5) = 6;
    dall.trlid(dall.trlid == 8) = 9;
    dall.trlid(dall.trlid == 4) = 7;
    d2ifc = dall;
    d2ifc(d2ifc.badRT ~= 0,:) = [];
    d2ifc(d2ifc.RT < 0.1,:) = [];
    
    % Use median response time per participant
    beh_stats_2ifc = grpstats(d2ifc,{'partid'},{'median', 'mean', 'range'},'datavars','RT');
    beh_stats_2ifc_2f = grpstats(d2ifc(d2ifc.trlid == 3,:),{'partid'},{'median', 'mean', 'range'},'datavars','RT');
    beh_stats_2ifc_2f1s = grpstats(d2ifc(d2ifc.trlid == 6,:),{'partid'},{'median', 'mean', 'range'},'datavars','RT');
    beh_stats_2ifc_1f2s = grpstats(d2ifc(d2ifc.trlid == 9,:),{'partid'},{'median', 'mean', 'range'},'datavars','RT');
    
    median_matrix = [beh_stats_2ifc_2f.median_RT, beh_stats_2ifc_2f1s.median_RT, ...
                     beh_stats_2ifc_1f2s.median_RT, beh_stats_2ifc.median_RT];
    
    % For now, exclude non-sensical / clearly bad fits
    % Having a threshold of larger than 0.5 means that even if the fit is
    % accepted, I cannot reasonably assume that I captured the PF properly,
    % because I am only measuring it's tail...
    bad_fit  = isnan(median_matrix) | repmat(~any(~isnan(peak_mat),2), [1, 4]);
    good_fit = ~bad_fit;
    
    
    % Correlate iAF with perceptual windows for all tasks (most importantly dfi)
    load(fullfile(data_dir, 'full_channel_vect.mat'));
    chanvect = full_channel_vect; clear full_channel_vect
    if ~exist(fullfile(fig_dir, fold_data), 'dir');
        mkdir(fullfile(fig_dir, fold_data));
    end
    percwin = ones(size(peak_mat))./peak_mat;
    col_vect = [[0.4 1 0.4]; [0.4 0.4 1]; [1 0.4 0.4]; [0.4 0.4 0.4]];
    [spearRho, pval, nobs, r, b1, b0] = deal(zeros(1, 4));
    fh = figure('color', 'w', 'position', [50 50 800, 790]);
    for ic = 1:4
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), median_matrix(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(percwin(good_fit(:,ic),1), median_matrix(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));

        subplot(2,2,ic)
        plot(percwin(good_fit(:,ic),1), median_matrix(good_fit(:,ic),ic), '.', 'color', col_vect(ic,:), 'markersize', 20); hold on
        line([1/12  1/7], [b0(1,ic)+b1(1,ic)*1/12 b0(1,ic)+b1(1,ic)*1/7], 'color', col_vect(ic,:)./2);
        xlim([1/12, 1/7]); ylim([min(min(median_matrix))*0.95 max(max(median_matrix))*1.05]);
        xticks = 0.07:0.01:0.14; set(gca, 'XTick', xticks)
        yticks = -0.1:0.1:1.8;   set(gca, 'YTick', yticks);
        addtext(sprintf('\nRho = %.2f\niAF = %.2f + %.2f*pWin\np-value = %.2f\nN = %.2f', spearRho(1,ic), b0(1,ic), b1(1,ic), pval(1,ic), nobs(1,ic)));
        if ic == 1
            title('2 Flashes vs 1 Flash');
        elseif ic == 2
            title('2 Flashes, 1 Sound vs 1 Flash, 1 Sound');
        elseif ic == 3
            title('1 Flash, 2 Sounds vs 2 Flashes, 2 Sounds');
        else
            title('All trials');
        end
        if ic > 2, xlabel('Occipital alpha cycle duration (s)'); end
        if ~isEven(ic), ylabel('Median response time (s)'); end
        plotspecs
    end
    suptitle(sprintf('iAF and response times'))
    [spearRho, pval, nobs]
    saveas(fh, fullfile(fig_dir, an_fold, fold_data, sprintf('yn_iAF_RTs_eyes_open.emf')))
    close all
    save(fullfile(fig_dir, an_fold, fold_data, 'yn_iAF_RTs_data'), 'spearRho', 'pval', 'r', 'b1', 'b0', 'nobs')
    
    
    
    
end % data subset (folder name)





% // eof

















