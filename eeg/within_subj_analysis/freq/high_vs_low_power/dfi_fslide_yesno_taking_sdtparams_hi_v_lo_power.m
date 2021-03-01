% Within subject frequency sliding analysis script.
% Yes-no task data.
% Focus on signal detection parameters (d' and bias). 
%
% Parent script(s): 
%   dfi_fslide_yesno_taking_sessions_into_account_sdtparams.m
%   dfi_tfrs_yesno_taking_sessions_into_account_sdtparams.m
%
% Children script(s): 
%   dfi_clusterPerm_dprime_freq_yesno_hi_lo_pow.m
%   dfi_clusterPerm_criterion_freq_yesno_hi_lo_pow.m
%
%
% Sibling script(s):
%   dfi_fslide_taking_sdtparams_hi_v_lo_power.m
%
%
% DETAILS
%
% The same analysis as dfi_fslide_yesno_taking_sessions_into_account_sdtparams.m, 
% but separately for high and low alpha power quantiles (median split). 
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


clc
close all
clear all


%% ------ USER INPUT -------

% Do we use high or low alpha power?

% 18.01.2021; 09:48h: Running hi_pow
% 18.01.2021; 10.00h: Running lo_pow


pow_cond = 'lo_pow';  % lo_pow, hi_pow


fprintf('\nSignal detection analysis of high versus low frequency prestimulus\n')
fprintf('focusing on %s trials\n\n', pow_cond);


%% 0.) --- SETUP ---


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
task     = 'yesno';

% Which channels do we want to investigate?
channels_of_interest = {'O2', 'PO4', 'PO8'};


%% Disentangle high from low power trials

% Since I have trial resolved data using the frequency sliding matrix
% fslide_mat, and the behavioural table btable, simply sub-select trials
% for high / low alpha power at the beginning and run the rest of the 
% script like before.


% Load power data
load(fullfile(data_dir, 'tfr_yesno_table_alpha_only.mat'), 'tfr_mat', 'btable');
btable_pow = btable; clear btable

% Load frequency sliding data
load(fullfile(data_dir, 'fslide_table_yesno.mat'), 'fslide_mat', 'btable');

% Make sure trials in tfr and fslide data match
assert(sum(btable_pow.trlid ~= btable.trlid) == 0)

% Regress out session effects etc. (before splitting into high / low pow)

% subject vector
partvect = unique(btable.partid);

% Channel of interest
choi = {'O2', 'PO4', 'PO8'};
[~,chid] = ismember(choi, channels_of_interest);

% Use linear model to account for variance due to time effects, because of
% session and trial number.

% Fit to average frequency in pre-stim window
% Get time axis
load(fullfile(data_dir, '701', 'yesno', 'inst_freq_yn', 'all_yn_sessions', ...
    'fslide_data_trls_all'));
toi = eeg_fslide.fslide_toi;
tid = toi > -0.602 & toi < -0.102;
btable.avgfslide = squeeze(nanmean(mean(fslide_mat(chid,tid,:),2),1));

% Make participant and sessid categorical variables
btable.sesscat = categorical(btable.sessid);
btable.partcat = categorical(btable.partid);
btable.runcat  = categorical(btable.run);

% Loop through subjects and fit separately, then remove effects from data

% Participant as fixed effect
lme2 = fitlm(btable,'avgfslide ~ sesscat * trl * block * partcat');
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
tfr_mat(trls_to_delete,:,:,:) = [];


% Get time and freq axes for power data
load(fullfile(data_dir, '701', 'yesno', 'tfr_zeropad', 'session_2', 'tfr_all_conds_all'));
toi = tfr_sub.time;
foi = tfr_sub.freq; clear tfr_sub

% Compute average alpha power in pre-stim window
choi = {'O2', 'PO4', 'PO8'};
[~,chid] = ismember(choi, channels_of_interest);
tid = toi > -0.602 & toi < -0.102;
alpha = foi > 5.9 & foi < 14.1;

% Add index column to btable
btable.idx = (1:length(btable))';

% Add average alpha power to btable
btable.avgpow_alpha = nanmean(squeeze(nanmean(nanmean(tfr_mat(:,chid,alpha,tid),3),4)),2);

% Divide into high and low power trial groups
[idx_hi, idx_lo] = deal([]);
for isubj = 1:length(partvect)
    
    % subj specific data
    partid = partvect(isubj);
    btable_subj = btable(btable.partid == partid,:);
    
    sessvect = unique(btable_subj.sessid);
    for isess = 1:length(sessvect)
        
        % sess specific data
        sessid = sessvect(isess);
        btable_sess = btable_subj(btable_subj.sessid == sessid,:);
        
        % median split
        
        % tfr_sess dimord:
        % trials        chans      freq        time
        % 1598           3          21          16
        [B,I] = sort(btable_sess.avgpow_alpha);
        
        % Keep track of trials with high / low alpha
        half = ceil(length(I)/2);
        idx_lo = [idx_lo; btable_sess.idx(I(1:half-1))];
        idx_hi = [idx_hi; btable_sess.idx(I(half:end))];
        
        %figure;
        %plot(1:length(idx_lo), btable.avgpow_alpha(idx_lo), '.b', 'markersize', 4); hold on
        %plot(length(idx_lo)+1:length(idx_lo)+length(idx_hi), btable.avgpow_alpha(idx_hi), '.r', 'markersize', 4)
        
    end % sess
    
end % subj

% Check if median split worked
fprintf('\n')
disp('SANITY CHECK - DID WE DIVIDE TRIALS INTO HIGH VS LOW POWER CORRECTLY?')
fprintf('\n')
fprintf('Avg lo pow = %f +- %f SEM\n', mean(btable.avgpow_alpha(idx_lo)), ...
    std(btable.avgpow_alpha(idx_lo)) / sqrt(20))
fprintf('Avg hi pow = %f +- %f SEM\n', mean(btable.avgpow_alpha(idx_hi)), ...
    std(btable.avgpow_alpha(idx_hi)) / sqrt(20))
fprintf('\n')



% Select only high or low power trials 
if strcmp(pow_cond, 'hi_pow')
    
   btable = btable(idx_hi, :);
   fslide_clean = fslide_clean(:, idx_hi);
   
elseif strcmp(pow_cond, 'lo_pow')
    
   btable = btable(idx_lo, :);
   fslide_clean = fslide_clean(:, idx_lo);  
   
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

disp('Computing dprime and bias for first and third tercile')
disp('of pre-stimulus alpha frequency.')
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
mkdir(fullfile(data_dir, 'sdt', 'freq_slide', pow_cond))
save(fullfile(data_dir, 'sdt', 'freq_slide', pow_cond, ...
    'sd_params_d_c_yesno_tcollapse.mat'), ...
    'dp_mat', 'dp_within_SE', 'c_mat', 'c_within_SE');




%% SDT approach over time (same 1F trials per SOA)

% Time windows over which to compute d' and C
toi = eeg_fslide.fslide_toi;
time_windows = -0.602:0.04:-0.078;

nt = length(time_windows);

% Divide data according to alpha power into nq quantiles
nq = 3; % 2 = median split

% Go through each SOA separately and compute SDT parameters
soas = unique(btable.soa);
[dp_cell, c_cell] = deal(cell(length(soas),nt-1));

disp('Computing dprime and bias for first and third tercile')
disp('of pre-stimulus alpha frequency for each time point.')
for isoa = 1:length(soas)
    
    fprintf('\nSOA %f2.0', soas(isoa))

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
        nq = 3; % 2 = median split
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

            end % quantile loop
        end % subject loop

        % Accumulate data over time windows:
        dp_cell{isoa,t} = dp_mat;
        c_cell{isoa,t} = c_mat;

    end % time windows

end % soas

% Is the false alarm rate always the same?
fa_mat(:,s,1,1);

% Accumulate continuous d' and criterion data in matrix
% levels of c_mat_cont and dp_mat_cont: soa x t x subj x cond x qntl 
dp_mat_cont = nan([length(soas),nt-1,size(dp_mat)]);
c_mat_cont  = nan([length(soas),nt-1,size(c_mat)]);
for isoa = 1:length(soas)
    for t = 1:nt-1
        dp_mat_cont(isoa,t,:,:,:) = dp_cell{isoa,t};
        c_mat_cont(isoa,t,:,:,:) = c_cell{isoa,t};
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
mkdir(fullfile(data_dir, 'sdt', 'freq_slide', pow_cond))
save(fullfile(data_dir, 'sdt', 'freq_slide', pow_cond, 'sd_params_d_c_yesno.mat'), ...
    'dp_mat_cont', 'dp_within_SE', 'c_mat_cont', 'c_within_SE', 'tvect', 'time_windows');


% enable warnings again
warning('on','all')


% eof

