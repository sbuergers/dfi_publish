% Within subject frequency sliding analysis script in source space.
% Yes-no threshold task data.
% Focus on signal detection parameters (d' and bias). 
%
% Parent script(s): 
%   dfi_source_fslide_prep.m
%
% Children script(s): 
%   dfi_source_clusterPerm_dprime_freq.m
%   dfi_source_clusterPerm_criterion_freq.m
%
% Sibling script(s):
%   dfi_source_fslide_yesno_taking_sessions_into_account_sdtparams.m
%
%
% DETAILS
%
% Sort data (instantaneous frequency) in ascending order for 
% a.) each time point, and b.) collapsed over time. Select the
% first and third tercile and compute d' and bias for each. 
%
% If the alpha temporal resolution hypothesis holds true, we would
% expect a higher temporal resolution as indicated by a larger dprime
% for the third tercile (where frequency is high). 
%
% Time window of interest: -0.6 to -0.1 seconds relative to stim-onset.
% Voxels of interest: See dfi_source_determine_func_roi.m
%
% Regress out nuisance factors such as session, trial etc...
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
an_fold = 'src_fslide_ynt';

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
save_dir = fullfile(data_dir, 'source_analysis', 'sdt', 'freq_slide');

% add fieldtrip folder to search path
try
    addpath(fullfile('fieldtrip-20160816'))
catch
    warning('Cannot find fieldtrip folder')
end

subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
    '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N        = 20;
task     = 'yn_threshold';

% Which channels do we want to investigate?
channels_of_interest = {'src1', 'src2', 'src3', 'src4', 'src5', 'src6', 'src7'};
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

        if isubj == 2 && isess == 2
            continue  % 702 session_6 should not be considered!
        end
        fprintf('Session %i\n', isess);

        sess_fold = sessnames{isess};
        load(fullfile(subj_dir, sess_fold, 'fslide_data_trls_all'));
        
        % Select channels
        [~, ch_indeces] = ismember(channels_of_interest, eeg_fslide.label);

        % Accumulate data in large matrix
        fslide_cell{isubj, isess} = eeg_fslide.fslide_med(ch_indeces,:,:);

        % Add additional session column, referring to ordinal number
        beh_stim.sessid = repmat(isess, [size(beh_stim,1),1]);

        % I also have to obtain the full behavioral data as a table
        btable = [btable; beh_stim];

        % Clean up a bit
        clear eeg_fslide beh_stim

    end % session combinations

    clear eeg_cell_stim beh_cell_stim

end % loop over participants


% Get time axis
load(fullfile(subj_dir, sess_fold, 'fslide_data_trls_all'));
toi = eeg_fslide.fslide_toi;

clear eeg_fslide

% Accumulate frequency sliding in matrix
fslide_mat = nan(length(channels_of_interest),size(toi,1),size(btable,1));
cmin = 1;
cmax = 0;
for isubj = 1:20
    for isess = 1:3
        if ~isempty(fslide_cell{isubj, isess})
            len = sum(btable.partid==str2double(subjvect{isubj}) & btable.sessid == isess);
            cmax = cmax + len;
            fslide_mat(:,:,cmin:cmax) = fslide_cell{isubj, isess};
            cmin = cmax + 1;
        end
    end
end

% Now I can use fslide_mat and btable, as they are both in table format
save(fullfile(data_dir, 'source_analysis', 'fslide_table.mat'), ...
    'fslide_mat', 'btable', 'channels_of_interest', 'toi', '-v7.3');






%% Look at temporal effects

load(fullfile(data_dir, 'source_analysis', 'fslide_table.mat'));

% Channel of interest
choi = channels_of_interest;
[~,chid] = ismember(choi, channels_of_interest);

% Use linear model to account for variance due to time effects, because of
% session and trial number.

% Fit to average frequency in pre-stim window
tid = toi > -0.602 & toi < -0.102;
btable.avgfslide = squeeze(nanmean(mean(fslide_mat(chid,tid,:),2),1));
% Make participant and sessid categorical variables
btable.sesscat = categorical(btable.sessid);
btable.partcat = categorical(btable.partid);
btable.runcat  = categorical(btable.run);

% Participant as fixed effect
lme2 = fitlm(btable,'avgfslide ~ sesscat * trl * block * partcat');
%lme2 = fitlm(btable,'avgfslide ~ sesscat + trl + block + partcat');
Resid2 = lme2.Residuals.Raw;
Pred2 = predict(lme2,btable);
nansum((btable.avgfslide - Resid2) - Pred2);


% Now apply correction to data and recompute frequency sliding effect
fslide_clean = squeeze(nanmean(fslide_mat(chid,:,:),1)) - repmat(Pred2', [size(toi,1), 1]);



%% Delete too fast trials
trls_to_delete = btable.RT < 0.1;
btable(trls_to_delete,:) = [];
fslide_mat(:,:,trls_to_delete) = [];
fslide_clean(:,trls_to_delete) = [];




% Get subject specific means per condition
% fslide_subj & fslide_GA rows (similar to 2,3,5,6,8,9):
% v1_1, v1_2, v2_1, v2_2, v1s1_1, v1s1_2, v2s1_1, v2s1_2, v1s2_1, v1s2_2, v2s2_1, v2s2_2
partvect = unique(btable.partid);
trlids = unique(btable.trlid);
[fslide_subj, fslide_subj_uncorr] = deal(nan(length(trlids)*2,20,size(toi,1)));
ntrls = nan(length(trlids)*2,20);
btable.see2 = btable.acc;
btable.see2(ismember(btable.trlid, [2,5,8])) = ~btable.see2(ismember(btable.trlid, [2,5,8]));
for isubj = 1:20
    for icond = 1:length(trlids)
        % see 1
        row_sel_see1 = btable.partid == partvect(isubj) & ...
            btable.trlid == trlids(icond) & ...
            ~btable.see2;
        ntrls((icond-1)*2+1,isubj) = sum(row_sel_see1);
        fslide_subj((icond-1)*2+1,isubj,:) = squeeze(mean(fslide_clean(:,row_sel_see1),2));
        fslide_subj_uncorr((icond-1)*2+1,isubj,:) = squeeze(nanmean(mean(fslide_mat(chid,:,row_sel_see1),3),1));
        % see 2
        row_sel_see2 = btable.partid == partvect(isubj) & ...
            btable.trlid == trlids(icond) & ...
            btable.see2;
        ntrls((icond-1)*2+2,isubj) = sum(row_sel_see2);
        fslide_subj((icond-1)*2+2,isubj,:) = squeeze(mean(fslide_clean(:,row_sel_see2),2));
        fslide_subj_uncorr((icond-1)*2+2,isubj,:) = squeeze(nanmean(mean(fslide_mat(chid,:,row_sel_see2),3),1));
    end
end

% get GA and standard error of the mean
[fslide_SE, fslide_GA, fslide_GA_uncorr, fslide_SE_uncorr] = deal(nan(length(trlids)*2,size(toi,1)));
for icond = 1:length(trlids)*2
    fslide_GA(icond,:) = nanmean(fslide_subj(icond,:,:),2);
    fslide_SE(icond,:) = nanstd(fslide_subj(icond,:,:)) ./ sqrt(20);
    fslide_GA_uncorr(icond,:) = nanmean(fslide_subj_uncorr(icond,:,:),2);
    fslide_SE_uncorr(icond,:) = nanstd(fslide_subj_uncorr(icond,:,:)) ./ sqrt(20);
end

% Plotting colours
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2

cnd_lbls = {'1F', '2F', '1F1S', '2F2S', '1F2S', '2F2S'};


% 1.) Plot fslide see versus not see per condition
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 1000, 600]);
ha = tight_subplot(2, 3,[0.1 0.075],[0.05],[0.05]);
for icond = 1:length(trlids)
    axes(ha(icond));
    tid = toi > -0.602 & toi < -0.102;
    tif = toi(tid);
    shadedErrorBar(tif', fslide_GA_uncorr((icond-1)*2+2,tid), fslide_SE((icond-1)*2+2,tid), {'-b','markerfacecolor', col2{icond}}, 1); hold on
    shadedErrorBar(tif', fslide_GA_uncorr((icond-1)*2+1,tid), fslide_SE((icond-1)*2+1,tid), {'-r','markerfacecolor', col1{icond}}, 1);
    % accentuate actual GA time courses
    plot(tif', fslide_GA_uncorr((icond-1)*2+2,tid),  'color', col2{icond}, 'linewidth', 1.5); hold on; % correct
    plot(tif', fslide_GA_uncorr((icond-1)*2+1,tid),  'color', col1{icond}, 'linewidth', 1.5);          % incorrect
    xlim([-0.6 -0.1])
    title(cnd_lbls{icond})
end


% 2.) Plot corrected fslide see versus not see per condition
fh2 = figure('color', [1 1 1], 'Position', [0, 0, 1000, 600]);
ha = tight_subplot(2, 3,[0.1 0.075],[0.05],[0.05]);
for icond = 1:length(trlids)
    axes(ha(icond));
    tid = toi > -0.602 & toi < -0.102;
    tif = toi(tid);
    shadedErrorBar(tif', fslide_GA((icond-1)*2+2,tid), fslide_SE((icond-1)*2+2,tid), {'-b','markerfacecolor', col2{icond}}, 1); hold on
    shadedErrorBar(tif', fslide_GA((icond-1)*2+1,tid), fslide_SE((icond-1)*2+1,tid), {'-r','markerfacecolor', col1{icond}}, 1);
    % accentuate actual GA time courses
    plot(tif', fslide_GA((icond-1)*2+2,tid),  'color', col2{icond}, 'linewidth', 1.5); hold on; % correct
    plot(tif', fslide_GA((icond-1)*2+1,tid),  'color', col1{icond}, 'linewidth', 1.5);          % incorrect
    xlim([-0.6 -0.1])
    title(cnd_lbls{icond})
end
close all




%% SDT approach

% Average over time (-0.6 to -0.1 s)
prestim = [-0.602, -0.078];
is_pre_stim = toi >= prestim(1) & toi <= prestim(2);
fs_clean_avg = nanmean(fslide_clean(is_pre_stim,:),1);

% Sort by frequency for each subject x condition
condvect = unique(btable.trlid);
[fs_clean_sorted, beh_sorted] = deal(cell(20,6));
for s = 1:N
    for c = 1:6
        % sort by frequency
        [fs_clean_sorted{s,c},sortidx] = sort(fs_clean_avg(1,btable.partid == partvect(s) & ...
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
        n = length(fs_clean_sorted{s,c});
        while mod(n,nq) > 0
            fs_clean_sorted{s,c}(floor(n/2)) = [];
            beh_sorted{s,c}(floor(n/2),:) = [];
            n = length(fs_clean_sorted{s,c});
        end
        % obtain quantiles
        for q = 1:nq
            qincr = n/nq;
            eeg_quantiles{s,c,q} = fs_clean_sorted{s,c}(1+(q-1)*qincr:q*qincr);
            beh_quantiles{s,c,q} = beh_sorted{s,c}(1+(q-1)*qincr:q*qincr,:);
        end
    end
end

% Create new behavioural table including quantile of alpha frequency
btable_quantiles = cell(nq,1);
for s = 1:N
    for c = 1:6
        for q = 1:nq
            beh_quantiles{s,c,q}.af = eeg_quantiles{s,c,q}';
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
        [~, sdm98] = dfi_SDM(beh,  9, 8, 0);
        
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
mkdir(save_dir);
save(fullfile(save_dir, 'sd_params_d_c_tcollapse.mat'), ...
    'dp_mat', 'dp_within_SE', 'c_mat', 'c_within_SE');





%% SDT approach over time

% Time windows over which to compute d' and C
time_windows = -0.602:0.04:-0.078;
nt = length(time_windows);

% Divide data according to alpha power into nq quantiles
nq = 3; % 2 = median split

[dp_cell, c_cell, freq_cell] = deal(cell(1,nt-1));
for t = 1:nt-1
    
    % Average over time window
    prestim = time_windows(t:t+1);
    is_pre_stim = toi >= prestim(1) & toi <= prestim(2);
    fs_clean_avg = nanmean(fslide_clean(is_pre_stim,:),1);
    
    % Sort by frequency for each subject x condition
    condvect = unique(btable.trlid);
    [fs_clean_sorted, beh_sorted] = deal(cell(20,6));
    for s = 1:N
        for c = 1:length(condvect)
            % sort by frequency
            [fs_clean_sorted{s,c},sortidx] = sort(fs_clean_avg(1,btable.partid == partvect(s) & ...
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
        for c = 1:length(condvect)
            % delete median until equal number of trials in each quantile
            n = length(fs_clean_sorted{s,c});
            while mod(n,nq) > 0
                fs_clean_sorted{s,c}(floor(n/2)) = [];
                beh_sorted{s,c}(floor(n/2),:) = [];
                n = length(fs_clean_sorted{s,c});
            end
            % obtain quantiles
            for q = 1:nq
                qincr = n/nq;
                eeg_quantiles{s,c,q} = fs_clean_sorted{s,c}(1+(q-1)*qincr:q*qincr);
                beh_quantiles{s,c,q} = beh_sorted{s,c}(1+(q-1)*qincr:q*qincr,:);
            end
        end
    end
    
    % Create new behavioural table including quantile of alpha frequency
    btable_quantiles = cell(nq,1);
    for s = 1:N
        for c = 1:length(condvect)
            for q = 1:nq
                beh_quantiles{s,c,q}.af = eeg_quantiles{s,c,q}';
                btable_quantiles{q} = [btable_quantiles{q}; beh_quantiles{s,c,q}];
            end
        end
    end
    
    % compute 2 equal variance gaussian SDT model for each quantile
    [dp_mat, c_mat] = deal(nan(N,3,nq)); % nsubj x cond x quantiles
    freq_mat = nan(N,6,nq);
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
            freq_mat(s,1,q) = nanmean(eeg_quantiles{s,1,q});
            freq_mat(s,2,q) = nanmean(eeg_quantiles{s,2,q});
            freq_mat(s,3,q) = nanmean(eeg_quantiles{s,3,q});
            freq_mat(s,4,q) = nanmean(eeg_quantiles{s,4,q});
            freq_mat(s,5,q) = nanmean(eeg_quantiles{s,5,q});
            freq_mat(s,6,q) = nanmean(eeg_quantiles{s,6,q});
            
        end % quantile loop
    end % subject loop
    
    % Accumulate data over time windows:
    dp_cell{t} = dp_mat;
    c_cell{t} = c_mat;
    freq_cell{t} = freq_mat;
    
end % time windows

% Accumulate continuous d' and criterion data in matrix
% levels of c_mat_cont and dp_mat_cont: t x subj x cond x qntl 
dp_mat_cont = nan([nt-1,size(dp_mat)]);
c_mat_cont  = nan([nt-1,size(c_mat)]);
freq_mat_cont = nan([nt-1,size(freq_mat)]);
for t = 1:nt-1
    dp_mat_cont(t,:,:,:) = dp_cell{t};
    c_mat_cont(t,:,:,:) = c_cell{t};
    freq_mat_cont(t,:,:,:) = freq_cell{t};
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
c_within_SE = nan(size(squeeze(dp_mat_cont(:,:,1,:)),1),3,3);
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
save(fullfile(save_dir, 'sd_params_d_c.mat'), ...
    'dp_mat_cont', 'dp_within_SE', 'c_mat_cont', 'c_within_SE', 'tvect', 'time_windows');







% plot alpha frequency quantiles
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 600, 800]);
ha = tight_subplot(2, 3,[0.05 0.03],[0.18],[0.18]);
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2
colavg = [0 0.6 0; 0 0 1; 1 0 0];
title_vect = {'1F', '2F', '1F1S', '2F1S', '1F2S', '2F2S'};

tvect = time_windows(2:nt)-diff(time_windows);
freq_within_SE = nan(size(squeeze(freq_mat_cont(:,:,1,:)),1),3,3);
for icond = 1:6
    % get Cousineau within subject SE for plotting:
    % Cancel out between subject variability by subtracting the subject
    % mean from each subject and then adding the grand mean
    data = squeeze(freq_mat_cont(:,:,icond,:));
    subj_mean = nanmean(data,3);
    grand_mean = nanmean(subj_mean,2);
    data_corr = data - repmat(subj_mean, [1,1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    freq_within_SE(:,icond,:) = nanstd(data_corr,0,2)./sqrt(size(data_corr,2));
    axes(ha(icond));
    plot(tvect,squeeze(nanmean(freq_mat_cont(:,:,icond,1),2)), 'color', 'r'); hold on;
    plot(tvect,squeeze(nanmean(freq_mat_cont(:,:,icond,3),2)), 'color', 'b'); hold on;
    shadedErrorBar(tvect, squeeze(nanmean(freq_mat_cont(:,:,icond,1),2)), freq_within_SE(:,icond,1), {'-r','markerfacecolor', 'r'}, 1); hold on
    shadedErrorBar(tvect, squeeze(nanmean(freq_mat_cont(:,:,icond,3),2)), freq_within_SE(:,icond,3), {'-b','markerfacecolor', 'b'}, 1); hold on
    title(title_vect{icond}); 
    ylim([-2 2]); xlim([-0.62 -0.09]);
    if icond == 1 || icond == 4, ylabel(sprintf('Frequency (Hz)')); end
    if icond == 1, legend('Quantile 1', 'Quantile 3'); end
    if icond > 3, xlabel('Time (s)'); end
    if icond ~= 1 && icond ~= 4, set(gca, 'yticklabel', []); end
end

% Compute effect size of difference
% effect size
data = squeeze(nanmean(freq_mat_cont(:,:,:,:),3));
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


% enable warnings again
warning('on','all')


% eof

