% Within subject frequency sliding analysis.
% Yes-no threshold task data.
% Compare pre-stim frequency for see1 or see2 flash outcomes.
%
% Parent script(s): 
%   dfi_fslide_taking_sessions_into_account_prep.m
%   dfi_fslide_taking_sessions_into_account_sdtparams.m
%
% Children script(s): 
%   ?
%
% Sibling script(s):
%   dfi_clusterPerm_dprime_freq_by_perc.m
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
clc
close all
clear all


% Name of directory to save things to
an_fold = 'sdt_params';

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
load_dir = fullfile(data_dir, 'freq_slide');
fig_dir  = 'dfi_experiment_figures';
save_dir = fullfile(fig_dir, 'Paper_figures', 'iAF', 'fslide', 'ynt', 'fslide');

% add fieldtrip folder to search path
try
    addpath(fullfile('fieldtrip-20160816'))
catch
    warning('Cannot find fieldtrip folder')
end

% useful variables
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N        = 20;
task     = 'yn_threshold';
foldercell = {'all_yn_sessions'};
ifold = 1;
fold_data  = foldercell{ifold};

foldername = 'Frequency_sliding';
if ~exist(fullfile(fig_dir, an_fold, foldername, fold_data, task), 'dir')
    mkdir(fullfile(fig_dir, an_fold, foldername, fold_data, task));
end



%% *** PREPARE DATA ***

% load data for d-prime
load(fullfile(load_dir, 'fslide_see1_v_see2.mat'));



%% *** STATISTICS ***

condvect = {'v2', 'fus', 'fis', 'av2'};

for icond = 1:length(condvect)

    % create fake fieldtrip structure
    tvect = tif;
    d_q1 = cell(20,1);
    d_q2 = cell(20,1);
    for isubj = 1:20
        if ~any(isnan(fslide_for_stats(:,isubj,icond,1)))
            d_q1{isubj}.time = tvect;
            d_q1{isubj}.fslide = fslide_for_stats(:,isubj,icond,1)';
        end
        if ~any(isnan(fslide_for_stats(:,isubj,icond,2)))
            d_q2{isubj}.time = tvect;
            d_q2{isubj}.fslide = fslide_for_stats(:,isubj,icond,2)';
        end
    end
    
    % check if we have proper data for each subject
    missing_ID = [];
    for i = 1:length(d_q1)
        if ~isfield(d_q1{i}, 'fslide') || ~isfield(d_q2{i}, 'fslide')
            missing_ID = [missing_ID, i];
        end
    end
    
    d_q1(missing_ID) = [];
    d_q2(missing_ID) = [];
    
    % Add labels (average over PO4, O2 and PO8):
    for isubj = 1:length(d_q1)
        d_q1{isubj}.label{1} = 'PO4_O2_PO8';
        d_q2{isubj}.label{1} = 'PO4_O2_PO8';    
        d_q1{isubj}.dimord = 'chan_time';
        d_q2{isubj}.dimord = 'chan_time';
    end
    
    % set up
    trltype  = 'all'; 
    param    = 'fslide'; % field to do the cluster permutation over ( e.g. powspctrm ) 
    chan     = 'PO4_O2_PO8';
    statistic= 'depsamplesT_bf';
    a        = 0.05;
    clustera = 0.05; % the two-side correction is done by cfg.correcttail (0.05 means 0.05 for two-tails)
    numperm  = 5000;
    method   = 'timelockanalysis';

    % show specs
    fprintf(['\n', '\n',  ...
             'Running cluster permutation with:', '\n', ...
            ['Method        = ', method], '\n', ...
            ['Trials        = ', trltype], '\n', ...
            ['Channels      = ', chan], '\n', ...
            ['Statistic     = ', statistic], '\n', ...
            ['Alpha level   = ', num2str(a)], '\n', ...
            ['Cluster alpha = ', num2str(clustera)], '\n', ...
            ['Permutations  = ', num2str(numperm), '\n', '\n']]);

    % % prepare neighbouring channel groups
    % cfg            = [];
    % cfg.feedback   = 'yes';
    % cfg.method     = 'triangulation';
    % cfg.layout     = 'EEG1010.lay';
    % neighbours     = ft_prepare_neighbours(cfg, fsd{1});
    neighbours = [];

    % Cluster permutation setup
    cfg = [];
    cfg.neighbours       = neighbours;
    cfg.parameter        = param;
    cfg.method           = 'montecarlo';
    cfg.statistic        = statistic; % With cfg.statistic = 'ft_statfun_actvsblT', 
                                      % we choose the so-called activation-versus-baseline 
                                      % T-statistic. This statistic compares the power in 
                                      % every sample (i.e., a (channel,frequency,time)-triplet) 
                                      % in the activation period with the corresponding 
                                      % time-averaged power (i.e., the average over the temporal 
                                      % dimension) in the baseline period. The comparison of the 
                                      % activation and the time-averaged baseline power is 
                                      % performed by means of a dependent samples T-statistic. 
    cfg.channel          = chan;  
    cfg.latency          = [-0.6  -0.1];
    cfg.avgovertime      = 'no';
    cfg.avgoverfreq      = 'no';
    cfg.avgoverchan      = 'no';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = clustera;
    cfg.clusterstatistic = 'maxsum';
    %cfg.clustercritval   = 
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.correcttail      = 'alpha';
    cfg.alpha            = a;       % if not correcttail = 'alpha' , adjust with bonferroni for two tails - 0.05/2
    cfg.numrandomization = numperm; % otherwise correcttail does this for us (relevant when using method montecarlo)

    % get design matrix
    subj = size(d_q1,1);
    design = zeros(2,2*subj);
    for i = 1:subj
      design(1,i) = i;
    end
    for i = 1:subj
      design(1,subj+i) = i;
    end
    design(2,1:subj)        = 1;
    design(2,subj+1:2*subj) = 2;

    cfg.design = design;
    cfg.uvar   = 1; % unit variable (i.e. Subject)
    cfg.ivar   = 2; % independent variable

    % run cluster permutation
    [stat] = ft_timelockstatistics(cfg, d_q1{:}, d_q2{:});

    bf = stat.bf;
    

    %% *** FIGURES ***

    % calculate the grand average for each condition
    cfg = [];
    cfg.channel   = 'PO4_O2_PO8';
    cfg.latency   = [-0.6  -0.1];
    cfg.parameter = 'fslide';
    GA_q1         = ft_timelockgrandaverage(cfg, d_q1{:});  
    GA_q2         = ft_timelockgrandaverage(cfg, d_q2{:});

    % Plot results (gray area denotes significant cluster)
    tif = GA_q1.time;
    fslide_GA = nan(size(tif,2),2);
    fslide_GA(:,1) = GA_q1.avg;
    fslide_GA(:,2) = GA_q2.avg;
    col1 = {[.2 .8 .2], [.2 .8 .2], [.5 .2 .9], [.2 .8 .2]}; % correct
    col2 = {[.5 .2 .9], [.5 .2 .9], [.2 .8 .2], [.5 .2 .9]}; % incorrect
    yl   = [min(min(fslide_GA(tif>-0.61&tif<0.01,:)))-0.1, max(max(fslide_GA(tif>-0.61&tif<0.01,:)))+0.1];
    xl   = [-0.6 -0.1];
    
    fh = figure('color', 'w', 'Position', [0 0 600 600]);
    ha = tight_subplot(2, 1,[0.01 0.01],[0.15],[0.15]);
    % Instantaneous frequency
    axes(ha(1));
    sem_f1_w = fslide_se(:,icond,1);
    sem_f2_w = fslide_se(:,icond,2);
    ylabel('Instantaneous frequency (Hz)');
    plot(tif, fslide_GA(:,2),  'color', col1{icond}, 'linewidth', 3); hold on; % correct
    plot(tif, fslide_GA(:,1),  'color', col2{icond}, 'linewidth', 3);          % incorrect   
    shadedErrorBar(tif, fslide_GA(:,2), sem_f2_w', {'-b','markerfacecolor', col2{icond}}, 0); hold on
    shadedErrorBar(tif, fslide_GA(:,1), sem_f1_w', {'-r','markerfacecolor', col1{icond}}, 0);
    gridxy(get(gca, 'xtick'), get(gca, 'ytick'), 'color', [0.8, 0.8, 0.8])
    xlim([-0.6 -0.1])
    % Bayes factors
    axes(ha(2));
    plot(tif, log10(bf),  'color', 'k', 'linewidth', 3);
    xlabel('Time from target onset (s)'); 
    ylabel('Bayes Factor (log10)');
    xlim([-0.6 -0.1])
    
    % Save data to remake figures later (for paper)
    mkdir(save_dir)
    save(fullfile(save_dir, sprintf('Bayes_factors_PO4_O2_PO8_%s.mat', condvect{icond})), ...
        'bf');
    
    % Save data to remake figures later (for paper)
    save(fullfile(save_dir, sprintf('figure_data_PO4_O2_PO8_%s.mat', condvect{icond})), ...
        'fslide_GA', 'tif', 'stat', 'sem_f1_w', 'sem_f2_w');
    
    % Save figures
    if ~exist(fullfile(fig_dir, an_fold, foldername, fold_data, task), 'dir')
        mkdir(fullfile(fig_dir, an_fold, foldername, fold_data, task))
    end
    close all
    
end


% eof
