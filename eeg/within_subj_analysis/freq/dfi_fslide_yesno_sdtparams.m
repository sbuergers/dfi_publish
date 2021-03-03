% Within subject frequency sliding analysis script.
% Yes-no task data.
% Focus on signal detection parameters (d' and bias). 
%
% Parent script(s): 
%   dfi_fslide_yesno_taking_sessions_into_account_prep.m
%
% Children script(s): 
%   dfi_clusterPerm_dprime_no_regress_freq_yesno.m
%   dfi_clusterPerm_criterion_no_regress_freq_yesno.m
%
% Sibling script(s):
%   dfi_fslide_sdtparams.m
%   dfi_fslide_taking_sessions_into_account_sdtparams.m
%   dfi_fslide_yesno_taking_sessions_into_account_sdtparams.m
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
% Channels of interest: PO4, PO8, O2
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
an_fold = 'inst_freq_yn';

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
save_dir = fullfile(data_dir, 'sdt', 'freq_slide', 'no_regress');

% add fieldtrip folder to search path
try
    addpath(fullfile('fieldtrip-20160816'))
catch
    warning('Cannot find fieldtrip folder')
end

subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
    '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N        = 20;
task     = 'yesno';

% Which channels do we want to investigate?
channels_of_interest = {'O1', 'O2', 'Oz', 'POz', 'PO4', 'PO8'};


% You can also directly get the data from D:\dfi_experiment_data\eeg_data\experiment
btable = [];

% loop over subjects
for isubj = 1:N
    
    fprintf('Subject number %i\n', isubj);
    
    subj_dir = fullfile(data_dir, subjvect{isubj}, task, an_fold);
    
    % Get session folder for this participant
    sessnames = {'all_yn_sessions'};
    
    % loop over sessions
    for isess = 1:length(sessnames)
        
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
fslide_mat = nan(6,size(toi,1),size(btable,1));
cmin = 1;
cmax = 0;
for isubj = 1:20
    if ~isempty(fslide_cell{isubj})
        len = sum(btable.partid==str2double(subjvect{isubj}));
        cmax = cmax + len;
        fslide_mat(:,:,cmin:cmax) = fslide_cell{isubj};
        cmin = cmax + 1;
    end
end

% Now I can use fslide_mat and btable, as they are both in table format
save(fullfile(data_dir, 'fslide_table_yesno.mat'), ...
    'fslide_mat', 'btable', 'channels_of_interest', 'toi', '-v7.3')



%% Look at temporal effects

load(fullfile(data_dir, 'fslide_table_yesno.mat'))

% subject vector
partvect = unique(btable.partid);

% Channel of interest
choi = {'O2', 'PO4', 'PO8'};
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


%% Delete too fast trials
trls_to_delete = btable.RT < 0.1;
btable(trls_to_delete,:) = [];
fslide_mat(:,:,trls_to_delete) = [];
Pred2(trls_to_delete) = [];


% Now do not apply correction to data
%fslide_clean = squeeze(nanmean(fslide_mat(chid,:,:),1)) - repmat(Pred2', [size(toi,1), 1]);
fslide_clean = squeeze(nanmean(fslide_mat(chid,:,:),1));

% Get subject specific means per condition
% fslide_subj & fslide_GA rows (similar to 2,3,5,6,8,9):
% v1_1, v1_2, v2_1, v2_2, v1s1_1, v1s1_2, v2s1_1, v2s1_2, v1s2_1, v1s2_2, v2s2_1, v2s2_2
[fslide_subj, fslide_subj_uncorr] = deal(nan(12,20,size(toi,1)));
trlids = unique(btable.trlid);
ntrls = nan(12,20);
btable.see2 = btable.acc;
btable.see2(ismember(btable.trlid, [2,5,8])) = ~btable.see2(ismember(btable.trlid, [2,5,8]));
for isubj = 1:20
    for icond = 1:6
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
[fslide_SE, fslide_GA, fslide_GA_uncorr, fslide_SE_uncorr] = deal(nan(12,size(toi,1)));
for icond = 1:12
    fslide_GA(icond,:) = nanmean(fslide_subj(icond,:,:),2);
    fslide_SE(icond,:) = nanstd(fslide_subj(icond,:,:)) ./ sqrt(20);
    fslide_GA_uncorr(icond,:) = nanmean(fslide_subj_uncorr(icond,:,:),2);
    fslide_SE_uncorr(icond,:) = nanstd(fslide_subj_uncorr(icond,:,:)) ./ sqrt(20);
end


%% SDT approach (same 1F trials per SOA)

% Time windows over which to compute d' and C
time_window = [-0.602, -0.078];
nt = 1;

% Divide data according to alpha power into nq quantiles
nq = 3; % 2 = median split

% Go through each SOA separately and compute SDT parameters
soas = unique(btable.soa);
[dp_cell, c_cell, freq_cell] = deal(cell(length(soas),nt));

for isoa = 1:length(soas)
    
    fprintf('\nSOA %f', soas(isoa))

    % Average over time window
    prestim = time_window;
    is_pre_stim = toi >= prestim(1) & toi <= prestim(2);
    fs_clean_avg = nanmean(fslide_clean(is_pre_stim,:),1);

    % Sort by frequency for each subject x condition
    condvect = unique(btable.trlid);
    [fslide_clean_sorted, beh_sorted] = deal(cell(20,6));
    for s = 1:N
        for c = 1:6
            % sort by frequency
            [fslide_clean_sorted{s,c},sortidx] = sort(fs_clean_avg(btable.partid == partvect(s) & ...
                btable.trlid == condvect(c) & ... 
                ((btable.soa == soas(isoa) & ismember(btable.trlid,[3,6,8,9])) | ...
                 ismember(btable.trlid, [2,5]))));
            % also sort behavioural responses accordingly
            beh_temp = btable(btable.partid == partvect(s) & btable.trlid == condvect(c) & ...
                ((btable.soa == soas(isoa) & ismember(btable.trlid,[3,6,8,9])) | ...
                 ismember(btable.trlid, [2,5])),:);
            beh_sorted{s,c} = beh_temp(sortidx,:);
        end
    end

    % do quantile split
    [eeg_quantiles, beh_quantiles] = deal(cell(20,6,nq));
    for s = 1:N
        for c = 1:6
            % delete median until equal number of trials in each quantile
            n = length(fslide_clean_sorted{s,c});
            while mod(n,nq) > 0
                fslide_clean_sorted{s,c}(floor(n/2)) = [];
                beh_sorted{s,c}(floor(n/2),:) = [];
                n = length(fslide_clean_sorted{s,c});
            end
            % obtain quantiles
            for q = 1:nq
                qincr = n/nq;
                eeg_quantiles{s,c,q} = fslide_clean_sorted{s,c}(1+(q-1)*qincr:q*qincr);
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
    freq_mat = nan(N,6,nq);
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

            fa_mat(isoa,s,1,q) = sdm23.pF_adj;
            fa_mat(isoa,s,2,q) = sdm56.pF_adj;
            fa_mat(isoa,s,3,q) = sdm98.pF_adj;

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
    dp_cell{isoa} = dp_mat;
    c_cell{isoa} = c_mat;
    freq_cell{isoa} = freq_mat;

end % soas

% Is the false alarm rate always the same?
fa_mat(:,s,1,1);

% Accumulate d' and criterion data in matrix
% levels of c_mat_cont and dp_mat_cont: soa x t x subj x cond x qntl 
dp_mat_soas = nan([length(soas),size(dp_mat)]);
c_mat_soas  = nan([length(soas),size(c_mat)]);
freq_mat_soas = nan([length(soas),size(freq_mat)]);
for isoa = 1:length(soas)
    dp_mat_soas(isoa,:,:,:) = dp_cell{isoa};
    c_mat_soas(isoa,:,:,:) = c_cell{isoa};
    freq_mat_soas(isoa,:,:,:) = freq_cell{isoa};
end

% Collapse over interm soas
interm_soas = 3:6;
dp_mat = squeeze(nanmean(dp_mat_soas(interm_soas,:,:,:), 1));
c_mat = squeeze(nanmean(c_mat_soas(interm_soas,:,:,:), 1));
freq_mat = squeeze(nanmean(freq_mat_soas(interm_soas,:,:,:), 1));


% sensitivity
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 600, 200]);
ha = tight_subplot(1, 3,[0.05 0.03],[0.18],[0.18]);
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2
colavg = [0 0.6 0; 0 0 1; 1 0 0];
title_vect = {'1F vs 2F', '1F1S vs 2F1S', '1F2S vs 2F2S'};
dp_within_SE = nan(3,3);
for icond = 1:3
    % get Cousineau within subject SE for plottinD:
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
    ylim([1 2.5]); xlim([0.5, 2.5]);
    if icond == 1, ylabel(sprintf('dPrime')); end
    xlabel('Quantile'); 
    if icond ~= 1, set(gca, 'yticklabel', []); end
end

% bias
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 600, 200]);
ha = tight_subplot(1, 3,[0.05 0.03],[0.18],[0.18]);
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2
colavg = [0 0.6 0; 0 0 1; 1 0 0];
title_vect = {'1F vs 2F', '1F1S vs 2F1S', '1F2S vs 2F2S'};
c_within_SE = nan(3,3);
for icond = 1:3
    % get Cousineau within subject SE for plottinD:
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
save(fullfile(save_dir, 'sd_params_d_c_yesno_tcollapse.mat'), ...
    'dp_mat', 'dp_within_SE', 'c_mat', 'c_within_SE');




%% SDT approach over time (same 1F trials per SOA)

% Time windows over which to compute d' and C
time_windows = -0.602:0.04:-0.078;

nt = length(time_windows);

% Divide data according to alpha power into nq quantiles
nq = 3; % 2 = median split

% Go through each SOA separately and compute SDT parameters
soas = unique(btable.soa);
[dp_cell, c_cell, freq_cell] = deal(cell(length(soas),nt-1));

for isoa = 1:length(soas)
    
    fprintf('\nSOA %f2', soas(isoa))

    for t = 1:nt-1

        fprintf('\nTime point %i', t)

        % Average over time window
        prestim = time_windows(t:t+1);
        is_pre_stim = toi >= prestim(1) & toi <= prestim(2);
        fs_clean_avg = nanmean(fslide_clean(is_pre_stim,:),1);

        % Sort by frequency for each subject x condition
        condvect = unique(btable.trlid);
        [fslide_clean_sorted, beh_sorted] = deal(cell(20,6));
        for s = 1:N
            for c = 1:6
                % sort by frequency
                [fslide_clean_sorted{s,c},sortidx] = sort(fs_clean_avg(btable.partid == partvect(s) & ...
                    btable.trlid == condvect(c) & ... 
                    ((btable.soa == soas(isoa) & ismember(btable.trlid,[3,6,8,9])) | ...
                     ismember(btable.trlid, [2,5]))));
                % also sort behavioural responses accordingly
                beh_temp = btable(btable.partid == partvect(s) & btable.trlid == condvect(c) & ...
                    ((btable.soa == soas(isoa) & ismember(btable.trlid,[3,6,8,9])) | ...
                     ismember(btable.trlid, [2,5])),:);
                beh_sorted{s,c} = beh_temp(sortidx,:);
            end
        end

        % do quantile split
        [eeg_quantiles, beh_quantiles] = deal(cell(20,6,nq));
        for s = 1:N
            for c = 1:6
                % delete median until equal number of trials in each quantile
                n = length(fslide_clean_sorted{s,c});
                while mod(n,nq) > 0
                    fslide_clean_sorted{s,c}(floor(n/2)) = [];
                    beh_sorted{s,c}(floor(n/2),:) = [];
                    n = length(fslide_clean_sorted{s,c});
                end
                % obtain quantiles
                for q = 1:nq
                    qincr = n/nq;
                    eeg_quantiles{s,c,q} = fslide_clean_sorted{s,c}(1+(q-1)*qincr:q*qincr);
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
                
                fa_mat(isoa,s,1,q) = sdm23.pF_adj;
                fa_mat(isoa,s,2,q) = sdm56.pF_adj;
                fa_mat(isoa,s,3,q) = sdm98.pF_adj;
                
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
        dp_cell{isoa,t} = dp_mat;
        c_cell{isoa,t} = c_mat;
        freq_cell{isoa,t} = freq_mat;

    end % time windows

end % soas

% Is the false alarm rate always the same?
fa_mat(:,s,1,1);

% Accumulate continuous d' and criterion data in matrix
% levels of c_mat_cont and dp_mat_cont: soa x t x subj x cond x qntl 
dp_mat_cont = nan([length(soas),nt-1,size(dp_mat)]);
c_mat_cont  = nan([length(soas),nt-1,size(c_mat)]);
freq_mat_cont = nan([length(soas),nt-1,size(freq_mat)]);
for isoa = 1:length(soas)
    for t = 1:nt-1
        dp_mat_cont(isoa,t,:,:,:) = dp_cell{isoa,t};
        c_mat_cont(isoa,t,:,:,:) = c_cell{isoa,t};
        freq_mat_cont(isoa,t,:,:,:) = freq_cell{isoa,t};
    end
end


% plot dP for 1st and 3rd alpha frequency quantile
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 1400, 600]);
ha = tight_subplot(3,length(soas), [0.01 0.01],[0.07],[0.04]);
for isoa = 1:length(soas)
    col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
    col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2
    colavg = [0 0.6 0; 0 0 1; 1 0 0];
    title_vect = {'1F vs 2F', '1F1S vs 2F1S', '1F2S vs 2F2S'};
    tvect = time_windows(2:nt)-diff(time_windows);
    dp_within_SE = nan(size(squeeze(dp_mat_cont(:,:,:,1,:)),2),3,3);
    for icond = 1:3
        % get Cousineau within subject SE for plottinD:
        % Cancel out between subject variability by subtracting the subject
        % mean from each subject and then adding the grand mean
        data = squeeze(dp_mat_cont(isoa,:,:,icond,:));
        subj_mean = nanmean(data,3);
        grand_mean = nanmean(subj_mean,2);
        data_corr = data - repmat(subj_mean, [1,1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
        dp_within_SE(:,icond,:) = nanstd(data_corr,0,2)./sqrt(size(data_corr,2));
        axes(ha(isoa+(icond-1)*length(soas)));
        plot(tvect,squeeze(nanmean(dp_mat_cont(isoa,:,:,icond,1),3)), 'color', col1{icond}); hold on;
        plot(tvect,squeeze(nanmean(dp_mat_cont(isoa,:,:,icond,3),3)), 'color', col2{icond}); hold on;
        shadedErrorBar(tvect, squeeze(nanmean(dp_mat_cont(isoa,:,:,icond,1),3)), dp_within_SE(:,icond,1), {'-r','markerfacecolor', col1{icond}}, 1); hold on
        shadedErrorBar(tvect, squeeze(nanmean(dp_mat_cont(isoa,:,:,icond,3),3)), dp_within_SE(:,icond,3), {'-b','markerfacecolor', col1{icond}}, 1); hold on
        if icond == 1, title(sprintf('%ims', round(1000*soas(isoa)))); end
        ylim([-1 4.5]); xlim([-0.62 -0.09]);
        if isoa == 1, ylabel(sprintf('%s',title_vect{icond})); end
        if icond == 1 && isoa == 1, legend('Quantile 1', 'Quantile 3'); end
        if icond == 3, xlabel('Time'); end
        if icond < 3, set(gca, 'xticklabel', []); end
        if isoa ~= 1, set(gca, 'yticklabel', []); end
    end
end

% For plotting we want to average standard error over the intermediate 4
% SOAs:
for icond = 1:3
    data = squeeze(nanmean(dp_mat_cont(3:6,:,:,icond,:),1));
    subj_mean = nanmean(data,3);
    grand_mean = nanmean(subj_mean,2);
    data_corr = data - repmat(subj_mean, [1,1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    dp_within_SE(:,icond,:) = nanstd(data_corr,0,2)./sqrt(size(data_corr,2));
end

% plot C for each alpha frequency quantile over time
% plot dP for 1st and 3rd alpha frequency quantile
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 1400, 600]);
ha = tight_subplot(3,length(soas), [0.01 0.01],[0.07],[0.04]);
for isoa = 1:length(soas)
    col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
    col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2
    colavg = [0 0.6 0; 0 0 1; 1 0 0];
    title_vect = {'1F vs 2F', '1F1S vs 2F1S', '1F2S vs 2F2S'};
    tvect = time_windows(2:nt)-diff(time_windows);
    c_within_SE = nan(size(squeeze(c_mat_cont(:,:,:,1,:)),2),3,3);
    for icond = 1:3
        % get Cousineau within subject SE for plottinD:
        % Cancel out between subject variability by subtracting the subject
        % mean from each subject and then adding the grand mean
        data = squeeze(c_mat_cont(isoa,:,:,icond,:));
        subj_mean = nanmean(data,3);
        grand_mean = nanmean(subj_mean,2);
        data_corr = data - repmat(subj_mean, [1,1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
        c_within_SE(:,icond,:) = nanstd(data_corr,0,2)./sqrt(size(data_corr,2));
        axes(ha(isoa+(icond-1)*length(soas)));
        plot(tvect,squeeze(nanmean(c_mat_cont(isoa,:,:,icond,1),3)), 'color', col1{icond}); hold on;
        plot(tvect,squeeze(nanmean(c_mat_cont(isoa,:,:,icond,3),3)), 'color', col2{icond}); hold on;
        shadedErrorBar(tvect, squeeze(nanmean(c_mat_cont(isoa,:,:,icond,1),3)), c_within_SE(:,icond,1), {'-r','markerfacecolor', col1{icond}}, 1); hold on
        shadedErrorBar(tvect, squeeze(nanmean(c_mat_cont(isoa,:,:,icond,3),3)), c_within_SE(:,icond,3), {'-b','markerfacecolor', col1{icond}}, 1); hold on
        if icond == 1, title(sprintf('%ims', round(1000*soas(isoa)))); end
        ylim([-1 3]); xlim([-0.62 -0.09]);
        if isoa == 1, ylabel(sprintf('%s',title_vect{icond})); end
        if icond == 1 && isoa == 1, legend('Quantile 1', 'Quantile 3', 'location', 'SouthEast'); end
        if icond == 3, xlabel('Time'); end
        if icond < 3, set(gca, 'xticklabel', []); end
        if isoa ~= 1, set(gca, 'yticklabel', []); end
    end
end

% For plotting we want to average standard error over the intermediate 4
% SOAs:
for icond = 1:3
    data = squeeze(nanmean(c_mat_cont(3:6,:,:,icond,:),1));
    subj_mean = nanmean(data,3);
    grand_mean = nanmean(subj_mean,2);
    data_corr = data - repmat(subj_mean, [1,1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    c_within_SE(:,icond,:) = nanstd(data_corr,0,2)./sqrt(size(data_corr,2));
end


% Save d-prime and criterion measures for statistics (cluster permutation
% test):
mkdir(save_dir)
save(fullfile(save_dir, 'sd_params_d_c_yesno.mat'), ...
    'dp_mat_cont', 'dp_within_SE', 'c_mat_cont', 'c_within_SE', 'tvect', 'time_windows');






% plot alpha frequency quantiles (for intermediate SOAs only)
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 1200, 1100]);
ha = tight_subplot(1, 6,[0.05 0.03],[0.18],[0.18]);
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]}; % see1
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]}; % see2
colavg = [0 0.6 0; 0 0 1; 1 0 0];
title_vect = {'1F', '2F', '1F1S', '2F1S', '1F2S', '2F2S'};

tvect = time_windows(2:nt)-diff(time_windows);
if size(freq_mat_cont,1) == 8
    freq_mat_cont = squeeze(nanmean(freq_mat_cont(3:6,:,:,:,:),1)); % average over interm. soas
end
freq_within_SE = nan(size(squeeze(freq_mat_cont(:,:,1,:)),1),3,3);
for icond = 1:6
    % get Cousineau within subject SE for plottinD:
    % Cancel out between subject variability by subtracting the subject
    % mean from each subject and then adding the grand mean
    data = squeeze(freq_mat_cont(:,:,icond,:));
    subj_mean = nanmean(data,3);
    grand_mean = nanmean(subj_mean,2);
    data_corr = data - repmat(subj_mean, [1,1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    freq_within_SE(:,icond,:) = nanstd(data_corr,0,2)./sqrt(size(data_corr,2));
    axes(ha(icond));
    plot(tvect,squeeze(nanmean(freq_mat_cont(:,:,icond,1),2)), 'color', 'r'); hold on;
    plot(tvect,squeeze(nanmean(freq_mat_cont(:,:,icond,2),2)), 'color', 'k'); hold on;
    plot(tvect,squeeze(nanmean(freq_mat_cont(:,:,icond,3),2)), 'color', 'b'); hold on;
    shadedErrorBar(tvect, squeeze(nanmean(freq_mat_cont(:,:,icond,1),2)), freq_within_SE(:,icond,1), {'-r','markerfacecolor', 'r'}, 1); hold on
    shadedErrorBar(tvect, squeeze(nanmean(freq_mat_cont(:,:,icond,2),2)), freq_within_SE(:,icond,2), {'-k','markerfacecolor', 'k'}, 1); hold on
    shadedErrorBar(tvect, squeeze(nanmean(freq_mat_cont(:,:,icond,3),2)), freq_within_SE(:,icond,3), {'-b','markerfacecolor', 'b'}, 1); hold on
    title(title_vect{icond}); 
    ylim([8.5 11]); xlim([-0.62 -0.09]);
    if icond == 1 || icond == 4, ylabel(sprintf('Frequency (Hz)')); end
    if icond == 1, legend('Quantile 1', 'Quantile 2', 'Quantile 3'); end
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




% plot alpha frequency quantiles for each subject, pooled over conditions
fh1 = figure('color', [1 1 1], 'Position', [0, 0, 1200, 1100]);
tvect = time_windows(2:nt)-diff(time_windows);
freq_mat_all_cond = squeeze(nanmean(freq_mat_cont,3));
for isubj = 1:N
    subplot(4,5,isubj)
    % get Cousineau within subject SE for plottinD:
    % Cancel out between subject variability by subtracting the subject
    % mean from each subject and then adding the grand mean
    data = squeeze(freq_mat_all_cond(:,isubj,:));
    subj_mean = nanmean(data,3);
    grand_mean = nanmean(subj_mean,2);
    data_corr = data - repmat(subj_mean, [1,1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    freq_within_SE(:,icond,:) = nanstd(data_corr,0,2)./sqrt(size(data_corr,2));
    plot(tvect, data(:,1), 'color', 'r'); hold on;
    plot(tvect, data(:,2), 'color', 'k'); hold on;
    plot(tvect, data(:,3), 'color', 'b'); hold on;
    title(sprintf('Subj %i',isubj)); 
    ylim([7 13]); xlim([-0.62 -0.09]);
    if icond == 1 || icond == 6 || icond == 11 || icond == 16, ylabel(sprintf('Frequency (Hz)')); end
    if icond == 1, legend('Quantile 1', 'Quantile 2', 'Quantile 3'); end
    if icond > 3, xlabel('Time (s)'); end
end
suptitle('Individual pre-stimulus frequency quantiles');


% enable warnings again
warning('on','all')


% eof

