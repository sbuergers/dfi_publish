% Within subject power analysis script.
% Yes-no threshold task data.
%
% Parent script(s): 
%   dfi_tfrs_zeropad_ynt_taking_sess_into_account_prep.m
%
% Children script(s): 
%   dfi_clusterPerm_dprime_no_regress_power.m
%   dfi_clusterPerm_criterion_no_regress_power.m
%
% Sibling script(s):
%   dfi_tfrs_yesno_sdtparams.m
%   dfi_tfrs_yesno_taking_sessions_into_account_sdtparams.m
%   dfi_tfrs_taking_sessions_into_account_sdtparams.m
%
%
% DETAILS
%
% Compute dprime and bias in pre-stim period over time for
% high and low power terciles.
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


%% 0.) --- SETUP ---
clc
close all
clear all

% disable warnings for speed sake:
warning('off','all')

% Name of directory to save things to
an_fold = 'tfr_zeropad';

% save figures?
save_figures = true;

% experiment script folder
try
    addpath(genpath('dfi'))
    dfi_startup
catch
    warning('Cannot find dfi folder')
end

% experiment data folder
data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
fig_dir  = 'dfi_experiment_figures';

% add fieldtrip folder to search path
try
    addpath(fullfile('fieldtrip-20160816'))
    addpath(fullfile('toolboxes', 'fieldtrip-20160816'))
catch
    warning('Cannot find fieldtrip folder')
end

subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
    '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N        = 20;
task     = 'yn_threshold';

% Which channels do we want to investigate?
channels_of_interest = {'O2', 'PO4', 'PO8'};
btable = [];

% loop over subjects
for isubj = 1:N
    
    fprintf('Subject number %i\n', isubj);
    
    subj_dir = fullfile(data_dir, subjvect{isubj}, task, an_fold);
    
    % Get session folder for this participant
    dirinfo = dir(subj_dir);
    dirnames = {dirinfo.name};
    Match = cellfun(@(x) strncmp(x, 'session_', 8), dirnames, 'UniformOutput', 0);
    sessnames = dirnames(cell2mat(Match));
    
    % loop over sessions
    for isess = 1:length(sessnames)
        
        fprintf('Session %i\n', isess);
        
        sess_fold = sessnames{isess};
        load(fullfile(subj_dir, sess_fold, 'tfr_all_conds_all'));
        
        % Remember time and frequency axes
        toi = tfr_sub.time;
        foi = tfr_sub.freq;
        
        % Select channels
        [~, ch_indeces] = ismember(channels_of_interest, tfr_sub.label);
        
        % Accumulate data in large matrix
        tfr_cell{isubj, isess} = tfr_sub.powspctrm(:,ch_indeces,:,:);
        
        % Add additional session column, referring to ordinal number
        beh_stim.sessid = repmat(isess, [size(beh_stim,1),1]);
        
        % I also have to obtain the full behavioral data as a table
        btable = [btable; beh_stim];
        
        % Clean up a bit
        clear tfr_sub beh_stim
        
    end % session combinations
    
    clear eeg_cell_stim beh_cell_stim
    
end % loop over participants


clear tfr_sub

% Accumulate frequency sliding in matrix
tfr_mat = nan(size(btable,1),3,size(foi,2),size(toi,2));
cmin = 1;
cmax = 0;
for isubj = 1:20
    for isess = 1:3
        if ~isempty(tfr_cell{isubj, isess})
            len = sum(btable.partid==str2double(subjvect{isubj}) & btable.sessid == isess);
            cmax = cmax + len;
            tfr_mat(cmin:cmax,:,:,:) = tfr_cell{isubj, isess};
            cmin = cmax + 1;
        end
    end
end

% Now I can use fslide_mat and btable, as they are both in table format
save(fullfile(data_dir, 'tfr_table_alpha_only.mat'), ...
    'tfr_mat', 'btable', 'channels_of_interest', 'toi', 'foi', '-v7.3');





%% Look at temporal effects

load(fullfile(data_dir, 'tfr_table_alpha_only.mat'), ...
    'tfr_mat', 'btable', 'channels_of_interest', 'toi', 'foi');


% subject vector
partvect = unique(btable.partid);

% Channel of interest
choi = {'O2', 'PO4', 'PO8'};
[~,chid] = ismember(choi, channels_of_interest);

% Use linear model to account for variance due to time effects, because of
% session and trial number.

% Fit to average frequency in pre-stim window
tid = toi > -0.602 & toi < -0.102;
alpha = foi > 5.9 & foi < 14.1;
btable.avgpow_alpha = nanmean(squeeze(nanmean(nanmean(tfr_mat(:,chid,alpha,tid),3),4)),2);


%% Delete too fast trials
trls_to_delete = btable.RT < 0.1;
btable(trls_to_delete,:) = [];
tfr_mat(trls_to_delete,:,:,:) = [];


% Make participant and sessid categorical variables
btable.sesscat = categorical(btable.sessid);
btable.partcat = categorical(btable.partid);
btable.runcat  = categorical(btable.run);

% Use log-transformed alpha power per trial
btable.log_alpha = log10(btable.avgpow_alpha);

% Participant as fixed effect
lme2 = fitlm(btable,'log_alpha ~ sesscat * trl * block * partcat');

Resid2 = lme2.Residuals.Raw;
Pred2 = predict(lme2,btable);
nansum((btable.avgpow_alpha - Resid2) - Pred2);


% Now apply correction to data and recompute power effect (in log space)
%tfr_clean = squeeze(nanmean(tfr_mat(:,chid,:,:),2)) - repmat(10.^(Pred2), [1, size(foi,2), size(toi,2)]);
tfr_clean = squeeze(nanmean(tfr_mat(:,chid,:,:),2));


% Get subject specific means per condition
% fslide_subj & fslide_GA rows (similar to 2,3,5,6,8,9):
% v1_1, v1_2, v2_1, v2_2, v1s1_1, v1s1_2, v2s1_1, v2s1_2, v1s2_1, v1s2_2, v2s2_1, v2s2_2
[tfr_subj, tfr_subj_uncorr] = deal(nan(12,20,size(foi,2),size(toi,2)));
trlids = unique(btable.trlid);
ntrls = nan(12,20);
btable.see2 = btable.acc;
btable.see2(ismember(btable.trlid, [2,5,8])) = ~btable.see2(ismember(btable.trlid, [2,5,8]));
for isubj = 1:20
    for icond = 1:6
        
        % Prepare trial indexing
        row_sel_see1 = btable.partid == partvect(isubj) & ...
            btable.trlid == trlids(icond) & ...
            ~btable.see2;
        
        row_sel_see2 = btable.partid == partvect(isubj) & ...
            btable.trlid == trlids(icond) & ...
            btable.see2;
        
        % Make sure to use the same number of trials
        n1 = sum(row_sel_see1);
        n2 = sum(row_sel_see2);
        ndiff = abs(n1-n2);
        if n1 > n2
            see1_indeces = find(row_sel_see1);
            row_sel_see1(see1_indeces(randsample([1:n1],ndiff))) = false;
        elseif n2 > n1
            see2_indeces = find(row_sel_see2);
            row_sel_see2(see2_indeces(randsample([1:n2],ndiff))) = false;
        end
        
        % see 1
        ntrls((icond-1)*2+1,isubj) = sum(row_sel_see1);
        tfr_subj((icond-1)*2+1,isubj,:,:) = squeeze(mean(tfr_clean(row_sel_see1,:,:),1));
        tfr_subj_uncorr((icond-1)*2+1,isubj,:,:) = squeeze(mean(nanmean(tfr_mat(row_sel_see1,chid,:,:),2),1));
        
        % see 2
        ntrls((icond-1)*2+2,isubj) = sum(row_sel_see2);
        tfr_subj((icond-1)*2+2,isubj,:,:) = squeeze(mean(tfr_clean(row_sel_see2,:,:),1));
        tfr_subj_uncorr((icond-1)*2+2,isubj,:,:) = squeeze(mean(nanmean(tfr_mat(row_sel_see2,chid,:,:),2),1));
    end
end

% get GA and standard error of the mean
[tfr_SE, tfr_GA, tfr_GA_uncorr, tfr_SE_uncorr] = deal(nan([size(tfr_subj,1), size(tfr_subj,3), size(tfr_subj,4)]));
for icond = 1:12
    tfr_GA(icond,:,:) = squeeze(nanmean(tfr_subj(icond,:,:,:),2));
    tfr_SE(icond,:,:) = squeeze(nanstd(tfr_subj(icond,:,:,:))) ./ sqrt(20);
    tfr_GA_uncorr(icond,:,:) = squeeze(nanmean(tfr_subj_uncorr(icond,:,:,:),2));
    tfr_SE_uncorr(icond,:,:) = squeeze(nanstd(tfr_subj_uncorr(icond,:,:,:))) ./ sqrt(20);
end
[tfr_GA_2vs1, tfr_SE_2vs1, tfr_GA_2vs1_uncorr, tfr_SE_2vs1_uncorr] = deal(nan([size(tfr_subj,1)/2, size(tfr_subj,3), size(tfr_subj,4)]));
for icond = 1:6
    tfr_GA_2vs1(icond,:,:) = squeeze(nanmean(tfr_subj((icond-1)*2+2,:,:,:) - tfr_subj((icond-1)*2+1,:,:,:),2));
    tfr_SE_2vs1(icond,:,:) = squeeze(nanstd(tfr_subj((icond-1)*2+2,:,:,:) - tfr_subj((icond-1)*2+1,:,:,:))) ./ sqrt(20);
    tfr_GA_2vs1_uncorr(icond,:,:) = squeeze(nanmean(tfr_subj_uncorr((icond-1)*2+2,:,:,:) - tfr_subj_uncorr((icond-1)*2+1,:,:,:),2));
    tfr_SE_2vs1_uncorr(icond,:,:) = squeeze(nanstd(tfr_subj_uncorr((icond-1)*2+2,:,:,:) - tfr_subj_uncorr((icond-1)*2+1,:,:,:))) ./ sqrt(20);
end

% Plotting colours
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2

cnd_lbls = {'1F', '2F', '1F1S', '2F1S', '1F2S', '2F2S'};
% 
% % Plot TFR See2-See1 per condition, TFR cleaned
% fh1 = figure('color', [1 1 1], 'Position', [0, 0, 1000, 600]);
% ha = tight_subplot(2, 3,[0.1 0.075],[0.05],[0.05]);
% for ci = 1:6
%     axes(ha(ci));
%     tfr_see2_min_see1 = squeeze(tfr_GA_2vs1(ci,:,:));
%     min_marg = min(min(tfr_see2_min_see1));
%     max_marg = max(max(tfr_see2_min_see1));
%     marg = max([abs(min_marg), abs(max_marg)]);
%     contourf(toi,foi,tfr_see2_min_see1,40,'linecolor','none')
%     set(gca,'clim',[-marg marg],'xlim',[-1.00 -0],'yscale','log','ytick',logspace(log10(min(foi)),log10(max(foi)),6),'yticklabel',round(logspace(log10(min(foi)),log10(max(foi)),6)*10)/10)
%     if ~(ci > 3), set(gca, 'xticklabel', []); end
%     if ci > 3, xlabel('Time'); end
%     title(sprintf('%s, see2-see1', cnd_lbls{ci}));
%     %colormap(haxby)
% end
% 
% % Plot TFR See2-See1 per condition, TFR not cleaned
% fh1 = figure('color', [1 1 1], 'Position', [0, 0, 1000, 600]);
% ha = tight_subplot(2, 3,[0.1 0.075],[0.05],[0.05]);
% for ci = 1:6
%     axes(ha(ci));
%     tfr_see2_min_see1 = squeeze(tfr_GA_2vs1_uncorr(ci,:,:));
%     min_marg = min(min(tfr_see2_min_see1));
%     max_marg = max(max(tfr_see2_min_see1));
%     marg = max([abs(min_marg), abs(max_marg)]);
%     contourf(toi,foi,tfr_see2_min_see1,40,'linecolor','none')
%     set(gca,'clim',[-marg marg],'xlim',[-1.00 -0],'yscale','log','ytick',logspace(log10(min(foi)),log10(max(foi)),6),'yticklabel',round(logspace(log10(min(foi)),log10(max(foi)),6)*10)/10)
%     if ~(ci > 3), set(gca, 'xticklabel', []); end
%     if ci > 3, xlabel('Time'); end
%     title(sprintf('%s, see2-see1', cnd_lbls{ci}));
%     %colormap(haxby)
% end



%% Individual peak frequency alpha power (iAP)
% Look at power for each trial closest to overal invididual trait alpha
% peak frequency

% Load iAF peak fits (using the toolbox written by Corcoran, 2017)
%load(fullfile(fig_dir, 'iAF_fits_corcoran', 'eyes_open_pkinfo.mat'));
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

% Get power only at individual alpha peak frequency
[iAP_clean, iAP] = deal(nan(size(tfr_clean,1), size(tfr_clean,3)));
tfr_chavg = squeeze(nanmean(tfr_mat,2));
for isubj = 1:20
    iAF = muPaf(isubj);
    [~, iAF_idx] = min(abs(foi-iAF));
    subj_idx = btable.partid == partvect(isubj);
    iAP_clean(subj_idx,:) = squeeze(tfr_clean(subj_idx,iAF_idx,:));
    iAP(subj_idx,:) = squeeze(tfr_chavg(subj_idx,iAF_idx,:));
end



%% SDT approach

% Average over time (-0.6 to -0.1 s)
prestim = [-0.602, -0.078];
is_pre_stim = toi >= prestim(1) & toi <= prestim(2);
iAP_clean_avg = nanmean(iAP_clean(:,is_pre_stim),2);

% Sort by frequency for each subject x condition
disp('Computing dprime and bias for first and third tercile')
disp('in pre-stimulus period.')
condvect = unique(btable.trlid);
[iAP_clean_sorted, beh_sorted] = deal(cell(20,6));
for s = 1:N
    for c = 1:6
        % sort by frequency
        [iAP_clean_sorted{s,c},sortidx] = sort(iAP_clean_avg(btable.partid == partvect(s) & ...
            btable.trlid == condvect(c)));
        % also sort behavioural responses accordingly
        beh_temp = btable(btable.partid == partvect(s) & btable.trlid == condvect(c),:);
        beh_sorted{s,c} = beh_temp(sortidx,:);
    end
end

% do quantile split
nq = 3; % 2 = median split
[eeg_quantiles, beh_quantiles] = deal(cell(20,6,nq));
for s = 1:N
    for c = 1:6
        % delete median until equal number of trials in each quantile
        n = length(iAP_clean_sorted{s,c});
        while mod(n,nq) > 0
            iAP_clean_sorted{s,c}(floor(n/2)) = [];
            beh_sorted{s,c}(floor(n/2),:) = [];
            n = length(iAP_clean_sorted{s,c});
        end
        % obtain quantiles
        for q = 1:nq
            qincr = n/nq;
            eeg_quantiles{s,c,q} = iAP_clean_sorted{s,c}(1+(q-1)*qincr:q*qincr);
            beh_quantiles{s,c,q} = beh_sorted{s,c}(1+(q-1)*qincr:q*qincr,:);
        end
    end
end

% Create new behavioural table including quantile of alpha frequency
btable_quantiles = cell(nq,1);
for s = 1:N
    for c = 1:6
        for q = 1:nq
            beh_quantiles{s,c,q}.ap = eeg_quantiles{s,c,q};
            btable_quantiles{q} = [btable_quantiles{q}; beh_quantiles{s,c,q}];
        end
    end
end

% compute 2 equal variance gaussian SDT model for each quantile
[dp_mat, c_mat] = deal(nan(N,3,nq)); % nsubj x cond x quantiles
for s = 1:N
    for q = 1:nq
        
        % 1F vs 2F
        beh = btable_quantiles{q}(btable_quantiles{q}.partid == partvect(s),:);
        
        % Type 1 signal detection parameters
        %S %N
        [~, sdm23] = dfi_SDM(beh,  3, 2, 0);
        [~, sdm56] = dfi_SDM(beh,  6, 5, 0);
        [~, sdm98] = dfi_SDM(beh,  8, 9, 0);
        
        dp_mat(s,1,q) = sdm23.dP_adj;
        dp_mat(s,2,q) = sdm56.dP_adj;
        dp_mat(s,3,q) = sdm98.dP_adj;
        
        c_mat(s,1,q) = sdm23.C_adj;
        c_mat(s,2,q) = sdm56.C_adj;
        c_mat(s,3,q) = sdm98.C_adj;
        
    end % quantile loop
end % subject loop


fh1 = figure('color', [1 1 1], 'Position', [0, 0, 600, 200]);
ha = tight_subplot(1, 3,[0.05 0.03],[0.18],[0.18]);
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2
colavg = [0 0.6 0; 0 0 1; 1 0 0];
title_vect = {'1F vs 2F', '1F1S vs 2F1S', '1F2S vs 2F2S'};
dp_within_SE = nan(3,3);
for icond = 1:3
    % get Cousineau within subject SE for plotting:
    % Cancel out between subject variability by subtracting the subject
    % mean from each subject and then adding the grand mean
    data = squeeze(dp_mat(:,icond,:));
    subj_mean = nanmean(data,2);
    grand_mean = nanmean(subj_mean,1);
    data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    dp_within_SE(icond,:) = nanstd(data_corr,0,1)./sqrt(size(data_corr,1));
    axes(ha(icond));
    plot(1,squeeze(nanmean(dp_mat(:,icond,1),1)), 'color', col1{icond}); hold on;
    plot(2,squeeze(nanmean(dp_mat(:,icond,3),1)), 'color', col2{icond}); hold on;
    errorbar(1, squeeze(nanmean(dp_mat(:,icond,1),1)), ...
         dp_within_SE(icond,1), 'color', colavg(icond,:), 'linewidth', 3); hold on
    errorbar(2, squeeze(nanmean(dp_mat(:,icond,3),1)), ...
         dp_within_SE(icond,3), 'color', colavg(icond,:), 'linewidth', 3); hold on
    title(title_vect{icond}); 
    ylim([0.5 2]); xlim([0.5, 2.5]);
    if icond == 1, ylabel(sprintf('dPrime')); end
    xlabel('Quantile'); 
    if icond ~= 1, set(gca, 'yticklabel', []); end
end


fh1 = figure('color', [1 1 1], 'Position', [0, 0, 600, 200]);
ha = tight_subplot(1, 3,[0.05 0.03],[0.18],[0.18]);
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2
colavg = [0 0.6 0; 0 0 1; 1 0 0];
title_vect = {'1F vs 2F', '1F1S vs 2F1S', '1F2S vs 2F2S'};
c_within_SE = nan(3,3);
for icond = 1:3
    % get Cousineau within subject SE for plotting:
    % Cancel out between subject variability by subtracting the subject
    % mean from each subject and then adding the grand mean
    data = squeeze(c_mat(:,icond,:));
    subj_mean = nanmean(data,2);
    grand_mean = nanmean(subj_mean,1);
    data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    c_within_SE(icond,:) = nanstd(data_corr,0,1)./sqrt(size(data_corr,1));
    axes(ha(icond));
    plot(1,squeeze(nanmean(c_mat(:,icond,1),1)), 'color', col1{icond}); hold on;
    plot(2,squeeze(nanmean(c_mat(:,icond,3),1)), 'color', col2{icond}); hold on;
    errorbar(1, squeeze(nanmean(c_mat(:,icond,1),1)), ...
         c_within_SE(icond,1), 'color', colavg(icond,:), 'linewidth', 3); hold on
    errorbar(2, squeeze(nanmean(c_mat(:,icond,3),1)), ...
         c_within_SE(icond,3), 'color', colavg(icond,:), 'linewidth', 3); hold on
    title(title_vect{icond}); 
    ylim([-1.2 1.2]);  xlim([0.5, 2.5]);
    if icond == 1, ylabel(sprintf('Bias')); end
    xlabel('Quantile'); 
    if icond ~= 1, set(gca, 'yticklabel', []); end
end


% Save d-prime and criterion measures
mkdir(fullfile(data_dir, 'sdt', 'power', 'no_regress'));
save(fullfile(data_dir, 'sdt', 'power', 'no_regress', 'sd_params_d_c_tcollapse.mat'), ...
    'dp_mat', 'dp_within_SE', 'c_mat', 'c_within_SE');



%% SDT approach over time

% Time windows over which to compute d' and C
time_windows = -0.602:0.04:-0.078;
nt = length(time_windows);

% Divide data according to alpha power into nq quantiles
nq = 3; % 2 = median split

disp('Computing dprime and bias for first and third tercile')
disp('time point by time point')
[dp_cell, c_cell, pow_cell] = deal(cell(1,nt-1));
for t = 1:nt-1
    
    fprintf('\nTime point %i', t)
    
    % Average over time window
    prestim = time_windows(t:t+1);
    is_pre_stim = toi >= prestim(1) & toi <= prestim(2);
    iAP_clean_avg = nanmean(iAP_clean(:,is_pre_stim),2);

    % Sort by frequency for each subject x condition
    condvect = unique(btable.trlid);
    [iAP_clean_sorted, beh_sorted] = deal(cell(20,6));
    for s = 1:N
        for c = 1:6
            % sort by frequency
            [iAP_clean_sorted{s,c},sortidx] = sort(iAP_clean_avg(btable.partid == partvect(s) & ...
                btable.trlid == condvect(c)));
            % also sort behavioural responses accordingly
            beh_temp = btable(btable.partid == partvect(s) & btable.trlid == condvect(c),:);
            beh_sorted{s,c} = beh_temp(sortidx,:);
        end
    end
    
    % do quantile split
    nq = 3; % 2 = median split
    [eeg_quantiles, beh_quantiles] = deal(cell(20,6,nq));
    for s = 1:N
        for c = 1:6
            % delete median until equal number of trials in each quantile
            n = length(iAP_clean_sorted{s,c});
            while mod(n,nq) > 0
                iAP_clean_sorted{s,c}(floor(n/2)) = [];
                beh_sorted{s,c}(floor(n/2),:) = [];
                n = length(iAP_clean_sorted{s,c});
            end
            % obtain quantiles
            for q = 1:nq
                qincr = n/nq;
                eeg_quantiles{s,c,q} = iAP_clean_sorted{s,c}(1+(q-1)*qincr:q*qincr);
                beh_quantiles{s,c,q} = beh_sorted{s,c}(1+(q-1)*qincr:q*qincr,:);
            end
        end
    end
    
    % Create new behavioural table including quantile of alpha frequency
    btable_quantiles = cell(nq,1);
    for s = 1:N
        for c = 1:6
            for q = 1:nq
                beh_quantiles{s,c,q}.af = eeg_quantiles{s,c,q};
                btable_quantiles{q} = [btable_quantiles{q}; beh_quantiles{s,c,q}];
            end
        end
    end
    
    % compute 2 equal variance gaussian SDT model for each quantile
    [dp_mat, c_mat] = deal(nan(N,3,nq)); % nsubj x cond x quantiles
    pow_mat = nan(N,6,nq);
    for s = 1:N
        for q = 1:nq
            
            % 1F vs 2F
            beh = btable_quantiles{q}(btable_quantiles{q}.partid == partvect(s),:);
            
            % Type 1 signal detection parameters
            %S %N
            [~, sdm23] = dfi_SDM(beh,  3, 2, 0);
            [~, sdm56] = dfi_SDM(beh,  6, 5, 0);
            [~, sdm98] = dfi_SDM(beh,  8, 9, 0);
            
            dp_mat(s,1,q) = sdm23.dP_adj;
            dp_mat(s,2,q) = sdm56.dP_adj;
            dp_mat(s,3,q) = sdm98.dP_adj;
            
            c_mat(s,1,q) = sdm23.C_adj;
            c_mat(s,2,q) = sdm56.C_adj;
            c_mat(s,3,q) = sdm98.C_adj;
            
            % keep track of frequencies as well (how variable is frequency
            % between quantiles anyway?)
            pow_mat(s,1,q) = nanmean(eeg_quantiles{s,1,q});
            pow_mat(s,2,q) = nanmean(eeg_quantiles{s,2,q});
            pow_mat(s,3,q) = nanmean(eeg_quantiles{s,3,q});
            pow_mat(s,4,q) = nanmean(eeg_quantiles{s,4,q});
            pow_mat(s,5,q) = nanmean(eeg_quantiles{s,5,q});
            pow_mat(s,6,q) = nanmean(eeg_quantiles{s,6,q});
            
        end % quantile loop
    end % subject loop
    
    % Accumulate data over time windows:
    dp_cell{t} = dp_mat;
    c_cell{t} = c_mat;
    pow_cell{t} = pow_mat;
    
end % time windows

% Accumulate continuous d' and criterion data in matrix
% levels of c_mat_cont and dp_mat_cont: t x subj x cond x qntl 
dp_mat_cont = nan([nt-1,size(dp_mat)]);
c_mat_cont  = nan([nt-1,size(c_mat)]);
pow_mat_cont = nan([nt-1,size(pow_mat)]);
for t = 1:nt-1
    dp_mat_cont(t,:,:,:) = dp_cell{t};
    c_mat_cont(t,:,:,:) = c_cell{t};
    pow_mat_cont(t,:,:,:) = pow_cell{t};
end


% plot dP for 1st and 3rd alpha frequency quantile
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 600, 200]);
ha = tight_subplot(1, 3,[0.05 0.03],[0.18],[0.18]);
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2
colavg = [0 0.6 0; 0 0 1; 1 0 0];
title_vect = {'1F vs 2F', '1F1S vs 2F1S', '1F2S vs 2F2S'};
tvect = time_windows(2:nt)-diff(time_windows);
dp_within_SE = nan(size(squeeze(dp_mat_cont(:,:,1,:)),1),3,3);
for icond = 1:3
    % get Cousineau within subject SE for plotting:
    % Cancel out between subject variability by subtracting the subject
    % mean from each subject and then adding the grand mean
    data = squeeze(dp_mat_cont(:,:,icond,:));
    subj_mean = nanmean(data,3);
    grand_mean = nanmean(subj_mean,2);
    data_corr = data - repmat(subj_mean, [1,1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    dp_within_SE(:,icond,:) = nanstd(data_corr,0,2)./sqrt(size(data_corr,2));
    axes(ha(icond));
    plot(tvect,squeeze(nanmean(dp_mat_cont(:,:,icond,1),2)), 'color', col1{icond}); hold on;
    plot(tvect,squeeze(nanmean(dp_mat_cont(:,:,icond,3),2)), 'color', col2{icond}); hold on;
    shadedErrorBar(tvect, squeeze(nanmean(dp_mat_cont(:,:,icond,1),2)), dp_within_SE(:,icond,1), {'-r','markerfacecolor', col1{icond}}, 1); hold on
    shadedErrorBar(tvect, squeeze(nanmean(dp_mat_cont(:,:,icond,3),2)), dp_within_SE(:,icond,3), {'-b','markerfacecolor', col1{icond}}, 1); hold on
    title(title_vect{icond}); 
    ylim([0.5 2]); xlim([-0.62 -0.09]);
    if icond == 1, ylabel(sprintf('dPrime')); end
    if icond == 1, legend('Quantile 1', 'Quantile 3'); end
    xlabel('Time'); 
    if icond ~= 1, set(gca, 'yticklabel', []); end
end


% plot C for each alpha frequency quantile over time
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 600, 200]);
ha = tight_subplot(1, 3,[0.05 0.03],[0.18],[0.18]);
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2
colavg = [0 0.6 0; 0 0 1; 1 0 0];
title_vect = {'1F vs 2F', '1F1S vs 2F1S', '1F2S vs 2F2S'};
tvect = time_windows(2:nt)-diff(time_windows);
c_within_SE = nan(size(data,1),3,3);
for icond = 1:3
    % get Cousineau within subject SE for plotting:
    % Cancel out between subject variability by subtracting the subject
    % mean from each subject and then adding the grand mean
    data = squeeze(c_mat_cont(:,:,icond,:));
    subj_mean = nanmean(data,3);
    grand_mean = nanmean(subj_mean,2);
    data_corr = data - repmat(subj_mean, [1,1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    c_within_SE(:,icond,:) = nanstd(data_corr,0,2)./sqrt(size(data_corr,2));
    axes(ha(icond));
    plot(tvect,squeeze(nanmean(c_mat_cont(:,:,icond,1),2)), 'color', col1{icond}); hold on;
    plot(tvect,squeeze(nanmean(c_mat_cont(:,:,icond,3),2)), 'color', col2{icond}); hold on;
    shadedErrorBar(tvect, squeeze(nanmean(c_mat_cont(:,:,icond,1),2)), c_within_SE(:,icond,1), {'-r','markerfacecolor', col1{icond}}, 1); hold on
    shadedErrorBar(tvect, squeeze(nanmean(c_mat_cont(:,:,icond,3),2)), c_within_SE(:,icond,3), {'-b','markerfacecolor', col1{icond}}, 1); hold on
    title(title_vect{icond}); 
    ylim([0.4 1.2]); xlim([-0.62 -0.09]);
    if icond == 1, ylabel(sprintf('Bias')); end
    if icond == 1, legend('Quantile 1', 'Quantile 3'); end
    xlabel('Time'); 
    if icond ~= 1, set(gca, 'yticklabel', []); end
end



% Save d-prime and criterion measures for statistics (cluster permutation
% test):
mkdir(fullfile(data_dir, 'sdt', 'power', 'no_regress'));
save(fullfile(data_dir, 'sdt', 'power', 'no_regress', 'sd_params_d_c.mat'), ...
    'dp_mat_cont', 'dp_within_SE', 'c_mat_cont', 'c_within_SE', 'tvect', 'time_windows');






% plot alpha frequency quantiles
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 1200, 1200]);
ha = tight_subplot(1, 6,[0.05 0.03],[0.18],[0.18]);
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2
colavg = [0 0.6 0; 0 0 1; 1 0 0];
title_vect = {'1F', '2F', '1F1S', '2F1S', '1F2S', '2F2S'};

tvect = time_windows(2:nt)-diff(time_windows);
pow_within_SE = nan(size(squeeze(pow_mat_cont(:,:,1,:)),1),3,3);
for icond = 1:6
    % get Cousineau within subject SE for plotting:
    % Cancel out between subject variability by subtracting the subject
    % mean from each subject and then adding the grand mean
    data = squeeze(pow_mat_cont(:,:,icond,:));
    subj_mean = nanmean(data,3);
    grand_mean = nanmean(subj_mean,2);
    data_corr = data - repmat(subj_mean, [1,1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    pow_within_SE(:,icond,:) = nanstd(data_corr,0,2)./sqrt(size(data_corr,2));
    axes(ha(icond));
    plot(tvect,squeeze(nanmean(pow_mat_cont(:,:,icond,1),2)), 'color', 'r'); hold on;
    plot(tvect,squeeze(nanmean(pow_mat_cont(:,:,icond,2),2)), 'color', 'k'); hold on;
    plot(tvect,squeeze(nanmean(pow_mat_cont(:,:,icond,3),2)), 'color', 'b'); hold on;
    shadedErrorBar(tvect, squeeze(nanmean(pow_mat_cont(:,:,icond,1),2)), pow_within_SE(:,icond,1), {'-r','markerfacecolor', 'r'}, 1); hold on
    shadedErrorBar(tvect, squeeze(nanmean(pow_mat_cont(:,:,icond,2),2)), pow_within_SE(:,icond,2), {'-k','markerfacecolor', 'k'}, 1); hold on
    shadedErrorBar(tvect, squeeze(nanmean(pow_mat_cont(:,:,icond,3),2)), pow_within_SE(:,icond,3), {'-b','markerfacecolor', 'b'}, 1); hold on
    title(title_vect{icond}); 
    ylim([-2 50]); xlim([-0.62 -0.09]);
    if icond == 1, ylabel(sprintf('Power')); end
    if icond == 1, legend('Quantile 1', 'Quantile 2', 'Quantile 3'); end
    xlabel('Time (s)');
    if icond ~= 1 && icond ~= 4, set(gca, 'yticklabel', []); end
end

% Compute effect size of difference
% effect size


data = squeeze(nanmean(pow_mat_cont(:,:,:,:),3));
cohens_d = nan(1,13);
for t = 1:13
    threshold_matrix = squeeze(data(t, :,:));
cohens_d(1,t) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,3) ) / ...
    sqrt(   nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,3)).^2 - ...
    2*corr(threshold_matrix(:,1), threshold_matrix(:,3)) * ...
    nanstd(threshold_matrix(:,1))*nanstd(threshold_matrix(:,3))   );
end
cohens_d
nanmean(cohens_d)





% plot alpha power quantiles for each subject, pooled over conditions
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 1200, 1100]);
tvect = time_windows(2:nt)-diff(time_windows);
pow_mat_all_cond = squeeze(nanmean(pow_mat_cont,3));
for isubj = 1:N
    subplot(4,5,isubj)
    % get Cousineau within subject SE for plotting:
    % Cancel out between subject variability by subtracting the subject
    % mean from each subject and then adding the grand mean
    data = squeeze(pow_mat_all_cond(:,isubj,:));
    plot(tvect, data(:,1), 'color', 'r'); hold on;
    plot(tvect, data(:,2), 'color', 'k'); hold on;
    plot(tvect, data(:,3), 'color', 'b'); hold on;
    title(sprintf('Subj %i',isubj)); 
    ylim([-1 200]); xlim([-0.62 -0.09]);
    if isubj == 1 || isubj == 6 || isubj == 11 || isubj == 16, ylabel(sprintf('Power')); end
    if isubj == 1, legend('Quantile 1', 'Quantile 2', 'Quantile 3'); end
    xlabel('Time (s)')
end
suptitle('Individual pre-stimulus power quantiles');




% plot alpha log-power quantiles for each subject, pooled over conditions
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 1200, 1100]);
tvect = time_windows(2:nt)-diff(time_windows);
pow_mat_all_cond = squeeze(nanmean(pow_mat_cont,3));
for isubj = 1:N
    subplot(4,5,isubj)
    % get Cousineau within subject SE for plotting:
    % Cancel out between subject variability by subtracting the subject
    % mean from each subject and then adding the grand mean
    data = squeeze(pow_mat_all_cond(:,isubj,:));
    plot(tvect, log10(data(:,1)), 'color', 'r'); hold on;
    plot(tvect, log10(data(:,2)), 'color', 'k'); hold on;
    plot(tvect, log10(data(:,3)), 'color', 'b'); hold on;
    title(sprintf('Subj %i',isubj)); 
    ylim([-1 3]); xlim([-0.62 -0.09]);
    if isubj == 1 || isubj == 6 || isubj == 11 || isubj == 16, ylabel(sprintf('log-Power')); end
    if isubj == 1, legend('Quantile 1', 'Quantile 2', 'Quantile 3'); end
    xlabel('Time (s)')
end
suptitle('Individual pre-stimulus log(power) quantiles');




% enable warnings again
warning('on','all')


% eof

