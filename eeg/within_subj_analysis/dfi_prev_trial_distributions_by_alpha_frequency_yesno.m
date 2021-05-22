% Control analysis: Look at trial distributions (conditions * SOAs)
% as a function of low or high pre-stimulus alpha frequency on the
% subsequent trial.
%
% Parent script(s): 
%   dfi_fslide_yesno_taking_sessions_into_account_prep.m
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
% Last modified May 2021


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

% Load fslide data
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


%% Delete too fast trials
trls_to_delete = btable.RT < 0.1;
btable(trls_to_delete,:) = [];


% For each participant and recording session, do a tercile split of alpha
% frequency and add a new factor to the tabular data (freq_tercile).
btable.freq_tercile = ones(size(btable,1), 1) * 2;
subjvect = unique(btable.partid);
for isubj = 1:length(subjvect)
    
    sessvect = unique(btable.sess(btable.partid == subjvect(isubj)));
    for isess = 1:length(sessvect)
        
        subj_sess_index = find(btable.partid==subjvect(isubj) & ...
                               btable.sess==sessvect(isess));
        lower_bound = quantile(btable.avgfslide(subj_sess_index), 0.33);
        upper_bound = quantile(btable.avgfslide(subj_sess_index), 0.67);
        btable.freq_tercile(...
            subj_sess_index(btable.avgfslide(subj_sess_index) < lower_bound) ...
        ) = 1;
        btable.freq_tercile(...
            subj_sess_index(btable.avgfslide(subj_sess_index) > upper_bound) ...
        ) = 3;
        
    end  % isess
    
end  % isubj


% For each participant, session, run and block, add the SOA and trlid of
% the previous trial
[btable.prev_trlid, btable.prev_soa] = deal(nan(size(btable,1), 1));
for isubj = 1:length(subjvect)
    
    subj_index = btable.partid == subjvect(isubj);
    sessvect = unique(btable.sess(subj_index));
    for isess = 1:length(sessvect)
        
        sess_index = btable.sess == sessvect(isess);
        runvect = unique(btable.run(subj_index & sess_index));
        for irun = 1:length(runvect)
            
            run_index = btable.run == runvect(irun);
            blockvect = unique(btable.run(subj_index & sess_index & run_index));
            for iblock = 1:length(blockvect)
                
                block_index = btable.block == blockvect(iblock);
                current_index = find(subj_index & sess_index & run_index & block_index);
                                
                prev_trlid = [nan; btable.trlid(current_index(2:end) - 1)];
                prev_soa = [nan; btable.soa(current_index(2:end) - 1)];
                
                btable.prev_trlid(current_index) = prev_trlid;
                btable.prev_soa(current_index) = prev_soa;
                
            end  % iblock
            
        end  % irun
        
    end  % isess
    
end  % isubj


% Distribution of trlids and soas of previous trials as a function of
% high or low alpha frequency on the current trial, for each subject
soavect = unique(btable.soa);
trlidvect = unique(btable.trlid);
[trlid_soa_counts_low, trlid_soa_counts_high] = deal( ...
    nan(length(subjvect), length(trlidvect), length(soavect)) ...
);
for isubj = 1:length(subjvect)
    
    bsubj = btable(btable.partid == subjvect(isubj), :);
    
    for itrlid = 1:length(trlidvect)
        
        for isoa = 1:length(soavect)
            
            trlid_soa_counts_low(isubj, itrlid, isoa) = sum( ...
                bsubj.freq_tercile == 1 & ...
                bsubj.prev_soa == soavect(isoa) & ...
                bsubj.prev_trlid == trlidvect(itrlid) ...
            );
            
            trlid_soa_counts_high(isubj, itrlid, isoa) = sum( ...
                bsubj.freq_tercile == 3 & ...
                bsubj.prev_soa == soavect(isoa) & ...
                bsubj.prev_trlid == trlidvect(itrlid) ...
            );
                
        end  % isoa
        
    end  % itrlid
    
end  % isubj

trlid_soa_count_ratio = trlid_soa_counts_high ./ (trlid_soa_counts_high + trlid_soa_counts_low);
trlid_soa_count_diff = trlid_soa_counts_high - trlid_soa_counts_low;

trlid_soa_count_diff_se = squeeze(std(trlid_soa_count_diff) ./ sqrt(20));
trlid_soa_count_diff_mean = squeeze(mean(trlid_soa_count_diff));


%% Figures

fh1 = figure('color', [1 1 1], 'Position', [0, 0, 1000, 705]);
ha = tight_subplot(6, 8, [0.02 0.02], [0.09], [0.09]);

i = 0;
condition_labels = {'1F', '2F', '1F1S', '2F1S', '1F2S', '2F2S'};
for itrlid = 1:length(trlidvect)
    for isoa = 1:length(soavect)
    
        i = i + 1;
    
        axes(ha(i));
        mean = trlid_soa_count_diff_mean(itrlid, isoa);
        stnd_err = trlid_soa_count_diff_se(itrlid, isoa);
        errorbar(0, mean, stnd_err, stnd_err, 'LineWidth', 1.5); hold on
        line([-10, 10], [0, 0])
        ylim([-3 3]);
        
        if i > 40
            xlabel(soavect(i-40));
        end

        if mod(i - 1, 8) == 0
            ylabel(condition_labels{itrlid});
        end
    end
end


% eof

