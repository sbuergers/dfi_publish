% TODO add scritp description and rename (looks at iAP not iAF)
% dfi_bb_TBW_consistency_over_tasks_and_conditions.m


%% *** SETUP ***
clear all

% Name of directory to save things to
an_fold = 'iAF_percwin_all_session_fit_betabinom';

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
figdir = 'D:\dfi_experiment_figures\PFs\beta_binom_weibull';

% add fieldtrip folder to search path
try
    addpath(fullfile('fieldtrip-20160816'))
catch
    warning('Cannot find fieldtrip folder')
end

% useful variables
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
filename = 'peaks_corcoran';
N = length(subjvect);

% Combine yn and ynt power estimates!
load(fullfile(fig_dir, 'iAF_fits_corcoran_zeropadded', 'eyes_open_pkinfo.mat'));
yn_pSpec = pSpec; clear pSpec
load(fullfile(fig_dir, 'iAF_fits_corcoran_zeropadded', 'yn_threshold', 'eyes_open_pkinfo.mat'));
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
[muPaf, muCog, muPow, muPow_osc, muPow_pk, muPow_pk_osc] = deal(nan(20,1));
for isubj = 1:20
    Nchans = [];
    for isess = 1:size(pSpec(isubj,:),2)
        if isstruct(pSpec(isubj,isess).chans)
            Nchans = [Nchans, size(pSpec(isubj,isess).chans,2)];
        end
    end
    [muPaf(isubj, :), muCog(isubj, :), muPow(isubj, :), muPow_osc(isubj, :), ...
        muPow_pk(isubj, :), muPow_pk_osc(isubj, :)] = meanIAF_sb([pSpec(isubj,:).sums], Nchans); 
end

% Which alpha power measure do you want to use?
alpha_power_variable = 'pow_abs'; % pow_osc, pkpow_abs, pkpow_osc

% Get peaks per condition for each subject and accumulate in matrix
foldercell = {'all_yn_and_ynt_sessions'};

% We might look at single sessions, or all of them combined (for now I only
% look at all of them combined, because I don't have psychometric function
% fits for single sessions...)
nfiles = 1;
for ifold = 1:nfiles
    
    % Prepare folders and peak matrices
    fold_data  = foldercell{ifold};
    if strcmp(alpha_power_variable, 'pow_abs')
        pow_mat = muPow;
    elseif strcmp(alpha_power_variable, 'pow_osc')
        pow_mat = muPow_osc;
    elseif strcmp(alpha_power_variable, 'pkpow_abs')
        pow_mat = muPow_pk;
    elseif strcmp(alpha_power_variable, 'pkpow_osc')
        pow_mat = muPow_pk_osc;
    end
    
    % make folders
    if ~exist(fullfile(fig_dir, an_fold, fold_data), 'dir')
        mkdir(fullfile(fig_dir, an_fold, fold_data));
    end

    %% *** YN-pooled THRESHOLD ***

    % Get psychometric function data
    % Here I already excluded bad fits. Plus separate function fits were
    % only substituted where the overall/joined fits failed!
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
    bad_fit  = isnan(ip) | repmat(~any(~isnan(pow_mat),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);
    
    % Plot correlations between perceptual window and iAF per channel and
    % condition
    percwin = pow_mat;
    col_vect = [[0.4 1 0.4]; [0.4 0.4 1]; [1 0.4 0.4]; [0.4 0.4 0.4]];
    [spearRho, pval, nobs, r, b1, b0] = deal(zeros(1, 4));
    fh = figure('color', 'w', 'position', [50 50 800, 790]);
    for ic = 1:3
        [pearsR(ic), pval_r(ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic));
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));

        subplot(2,2,ic)
        plot(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), '.', 'color', col_vect(ic,:), 'markersize', 20); hold on
        xl = xlim; yl = ylim;
        line([xl(1)  xl(2)], [b0(1,ic)+b1(1,ic)*xl(1) b0(1,ic)+b1(1,ic)*xl(2)], 'color', col_vect(ic,:)./2);
        addtext(sprintf('\nr^2 = %.2f\niAF = %.2f + %.2f*pWin\np-value = %.2f\nN = %.2f', spearRho(1,ic), b0(1,ic), b1(1,ic), pval(1,ic), nobs(1,ic)));
        if ic == 1
            title('2 Flashes');
        elseif ic == 2
            title('2 Flashes, 1 Sound');
        elseif ic == 3
            title('1 Flash, 2 Sounds');
        else
            title('2 Flashes, 2 Sounds');
        end
        if ic > 2, xlabel('Power (a.u.)'); end
        if ~isEven(ic), ylabel('Inflection point (s)'); end
        plotspecs
    end

    if ~exist(fullfile(fig_dir, an_fold, fold_data), 'dir');
        mkdir(fullfile(fig_dir, an_fold, fold_data));
    end
    export_fig(fh, fullfile(fig_dir, an_fold, fold_data, sprintf('yn_pooled_iAPow_perWin_corr_eyes_open_%s', alpha_power_variable)), '-tiff', '-m1.5')
    close all
    save(fullfile(fig_dir, an_fold, fold_data, sprintf('iAPow_perWin_data_yn_pooled_eyes_open_%s', alpha_power_variable)), 'spearRho', 'pval', 'pval_r', 'r', 'b1', 'b0', 'nobs')
    
    
    
    %% *** YN-pooled SLOPE ***

    % Get psychometric function data
    % Here I already excluded bad fits. Plus separate function fits were
    % only substituted where the overall/joined fits failed!
    folder = 'yn_pooled';
    load(fullfile(figdir, folder, 'data_betabinom_in_pffit_format_yesno.mat'));
    pffit_all = pffit_yn;
    
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
    percwin = pow_mat;
    col_vect = [[0.4 1 0.4]; [0.4 0.4 1]; [1 0.4 0.4]; [0.4 0.4 0.4]];
    [spearRho, pval, nobs, r, b1, b0] = deal(zeros(1, 4));
    fh = figure('color', 'w', 'position', [50 50 800, 790]);
    for ic = 1:3
        [pearsR(ic), pval_r(ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic));
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));

        subplot(2,2,ic)
        plot(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), '.', 'color', col_vect(ic,:), 'markersize', 20); hold on
        xl = xlim; yl = ylim;
        line([xl(1)  xl(2)], [b0(1,ic)+b1(1,ic)*xl(1) b0(1,ic)+b1(1,ic)*xl(2)], 'color', col_vect(ic,:)./2);
        addtext(sprintf('\nr^2 = %.2f\niAF = %.2f + %.2f*pWin\np-value = %.2f\nN = %.2f', spearRho(1,ic), b0(1,ic), b1(1,ic), pval(1,ic), nobs(1,ic)));
        if ic == 1
            title('2 Flashes');
        elseif ic == 2
            title('2 Flashes, 1 Sound');
        elseif ic == 3
            title('1 Flash, 2 Sounds');
        else
            title('2 Flashes, 2 Sounds');
        end
        if ic > 2, xlabel('Power (a.u.)'); end
        if ~isEven(ic), ylabel('log-slope'); end
        plotspecs
    end

    if ~exist(fullfile(fig_dir, an_fold, fold_data), 'dir');
        mkdir(fullfile(fig_dir, an_fold, fold_data));
    end
    export_fig(fh, fullfile(fig_dir, an_fold, fold_data, sprintf('yn_pooled_iAPow_slope_corr_eyes_open_%s', alpha_power_variable)), '-tiff', '-m1.5')
    close all
    save(fullfile(fig_dir, an_fold, fold_data, sprintf('iAPow_slope_data_yn_pooled_eyes_open_%s', alpha_power_variable)), 'spearRho', 'pval_r', 'pval', 'r', 'b1', 'b0', 'nobs')
    
    
    
    
    
    
    %% *** 2IFC THRESHOLD ***

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
    bad_fit  = isnan(ip) | repmat(~any(~isnan(pow_mat),2), [1, 3]);
    good_fit = ~bad_fit & ~(ip > 0.5);
    
    % Plot correlations between perceptual window and iAF per channel and
    % condition
    percwin = pow_mat;
    col_vect = [[0.4 1 0.4]; [0.4 0.4 1]; [1 0.4 0.4]; [0.4 0.4 0.4]];
    [spearRho, pval, nobs, r, b1, b0] = deal(zeros(1, 4));
    fh = figure('color', 'w', 'position', [50 50 800, 790]);
    for ic = 1:3
        [pearsR(ic), pval_r(ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic));
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));

        subplot(2,2,ic)
        plot(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic), '.', 'color', col_vect(ic,:), 'markersize', 20); hold on
        xl = xlim; yl = ylim;
        line([xl(1)  xl(2)], [b0(1,ic)+b1(1,ic)*xl(1) b0(1,ic)+b1(1,ic)*xl(2)], 'color', col_vect(ic,:)./2);
        addtext(sprintf('\nr^2 = %.2f\niAF = %.2f + %.2f*pWin\np-value = %.2f\nN = %.2f', spearRho(1,ic), b0(1,ic), b1(1,ic), pval(1,ic), nobs(1,ic)));
        if ic == 1
            title('2 Flashes');
        elseif ic == 2
            title('2 Flashes, 1 Sound');
        elseif ic == 3
            title('1 Flash, 2 Sounds');
        else
            title('2 Flashes, 2 Sounds');
        end
        if ic > 2, xlabel('Power (a.u.)'); end
        if ~isEven(ic), ylabel('Inflection point (s)'); end
        plotspecs
    end

    if ~exist(fullfile(fig_dir, an_fold, fold_data), 'dir');
        mkdir(fullfile(fig_dir, an_fold, fold_data));
    end
    export_fig(fh, fullfile(fig_dir, an_fold, fold_data, sprintf('twoifc_iAPow_perWin_corr_eyes_open_%s', alpha_power_variable)), '-tiff', '-m1.5')
    close all
    save(fullfile(fig_dir, an_fold, fold_data, sprintf('iAPow_perWin_data_twoifc_eyes_open_%s', alpha_power_variable)), 'spearRho', 'pval', 'pval_r', 'r', 'b1', 'b0', 'nobs')
    

    
    
    %% *** 2IFC SLOPE ***

    % Get psychometric function data
    % Here I already excluded bad fits. Plus separate function fits were
    % only substituted where the overall/joined fits failed!
    folder = '2ifc';
    load(fullfile(figdir, folder, 'data_betabinom_in_pffit_format_2IFC.mat'));
    pffit_all = pffit_2ifc;
    
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
    percwin = pow_mat;
    col_vect = [[0.4 1 0.4]; [0.4 0.4 1]; [1 0.4 0.4]; [0.4 0.4 0.4]];
    [spearRho, pval, nobs, r, b1, b0] = deal(zeros(1, 4));
    fh = figure('color', 'w', 'position', [50 50 800, 790]);
    for ic = 1:3
        [pearsR(ic), pval_r(ic)] = corr(percwin(good_fit(:,ic),1), ip(good_fit(:,ic),ic));
        [spearRho(1,ic), pval(1,ic)] = corr(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), 'type', 'Spearman', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), 'one');
        nobs(1,ic) = sum(good_fit(:,ic));

        subplot(2,2,ic)
        plot(percwin(good_fit(:,ic),1), slope(good_fit(:,ic),ic), '.', 'color', col_vect(ic,:), 'markersize', 20); hold on
        xl = xlim; yl = ylim;
        line([xl(1)  xl(2)], [b0(1,ic)+b1(1,ic)*xl(1) b0(1,ic)+b1(1,ic)*xl(2)], 'color', col_vect(ic,:)./2);
        addtext(sprintf('\nr^2 = %.2f\niAF = %.2f + %.2f*pWin\np-value = %.2f\nN = %.2f', spearRho(1,ic), b0(1,ic), b1(1,ic), pval(1,ic), nobs(1,ic)));
        if ic == 1
            title('2 Flashes');
        elseif ic == 2
            title('2 Flashes, 1 Sound');
        elseif ic == 3
            title('1 Flash, 2 Sounds');
        else
            title('2 Flashes, 2 Sounds');
        end
        if ic > 2, xlabel('Power (a.u.)'); end
        if ~isEven(ic), ylabel('log-slope'); end
        plotspecs
    end

    if ~exist(fullfile(fig_dir, an_fold, fold_data), 'dir');
        mkdir(fullfile(fig_dir, an_fold, fold_data));
    end
    export_fig(fh, fullfile(fig_dir, an_fold, fold_data, sprintf('twoifc_iAPow_slope_corr_eyes_open_%s', alpha_power_variable)), '-tiff', '-m1.5')
    close all
    save(fullfile(fig_dir, an_fold, fold_data, sprintf('iAPow_slope_data_twoifc_eyes_open_%s', alpha_power_variable)), 'spearRho', 'pval', 'pval_r', 'r', 'b1', 'b0', 'nobs')
    
    
    
    
    
    
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
    bad_fit  = isnan(pow_mat);
    good_fit = repmat(~(bad_fit), [1,4]);
    
    ip = soamat;
    
    % Correlate iAF with perceptual windows for all tasks (most importantly dfi)
    percwin = pow_mat;
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
        xl = xlim; yl = ylim;
        line([xl(1)  xl(2)], [b0(ic)+b1(ic)*xl(1) b0(ic)+b1(ic)*xl(2)], 'color', col_vect(ic,:)./2);
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
        if ic > 2, xlabel('Power (a.u.)'); end
        if ~isEven(ic), ylabel('Inflection point (s)'); end
        plotspecs
    end
    [pearsR, pval, nobs]
    foldername = 'ynt_soa';
    if ~exist(fullfile(fig_dir, an_fold, fold_data, foldername), 'dir');
        mkdir(fullfile(fig_dir, an_fold, fold_data, foldername));
    end
    export_fig(fh, fullfile(fig_dir, an_fold, fold_data, foldername, sprintf('iAPow_vs_ynt_soa_%s', alpha_power_variable)), '-tiff', '-m1.5')
    save(fullfile(fig_dir, an_fold, fold_data, foldername, sprintf('yn_iAPow_perWin_data_%s', alpha_power_variable)), 'spearRho', 'pval', 'pval_r', 'r', 'b1', 'b0', 'nobs')
    close all
    
    
end % data subset (folder name)


% eof

