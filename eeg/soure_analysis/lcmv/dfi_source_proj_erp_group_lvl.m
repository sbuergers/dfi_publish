% Group level contrasts of source power
%
% Parent script(s): 
%   dfi_source_proj_erp.m
%
% Children script(s): 
%   dfi_source_determine_func_roi.m
%
% Sibling script(s):
%   None
%
%
% DETAILS
%
% Compute group level contrasts of source power of the
% form (VAR(post)-VAR(pre)) ./ VAR(pre).
% (see dfi_source_proj_erp.m for details)
%
% NOTES
% For this source analysis we use a newer fieldtrip version than for
% the rest of the analyses (quite a bit has changed in fieldtrip 
% regarding sourceanalysis over the years so this simply seems 
% least likely to lead to problems).
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


% Set paths and load data
restoredefaultpath; clc; close all; clear all

% Add paths
try
    addpath(genpath('dfi'))
    dfi_startup
catch
    warning('Cannot find dfi folder')
end
ft_path = fullfile('fieldtrip20200906', 'fieldtrip-master'); 
addpath(ft_path);
ft_defaults

% Get data
data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
dirD = fullfile(data_dir, 'source_analysis');

% Useful variables
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
    '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N = length(subjvect);

% Get anatomical image (Colin27)
mri = ft_read_mri(fullfile(ft_path, 'template', 'headmodel', 'standard_mri.mat'));

% Get number of leadfield voxels in brain
load(fullfile(dirD, '701', 'sess2', ...
              'inv_sols_lcmv_c27_100to300ms_600to100msBl_1F2F.mat'), ...
     'leadfield');
n1 = sum(leadfield.inside);


%% Main loop
[effects_lcmv_subj, erp_subj_pre, erp_subj_post] = deal(cell(N, 1));
[ga_pre_lcmv, ga_post_lcvm, ga_lcmv, ga_pre_elor, ga_post_elor_simul, ...
    contrast2, contrast2_norm_sub, contrast2_norm_div, ...
    contrast2_mean, contrast2_mean_norm_sub, contrast2_mean_norm_div] = ...
    deal(cell(20,1));

for isubj = 1:N
    
    % Get session data
    partid = subjvect{isubj};
    
    fprintf('\n-----------------------------------')
    fprintf('\nParticipant %s\n', partid)
    fprintf('\n         Available sessions \n')
    ls(fullfile(dirD, sprintf('%s\\sess*',partid)))
    fprintf('-----------------------------------\n')
    
    sessions = ls(fullfile(dirD, sprintf('%s\\sess*',partid)));
    
    [effects_lcmv, erp_pre_sess, erp_post_sess, ...
        contrast_sess, contrast_sess_norm_sub, contrast_sess_norm_div, ...
        contrast_sess_mean, contrast_sess_mean_norm_sub, contrast_sess_mean_norm_div, ...
        lcmv_pre_cell, lcmv_post_cell, erp_post_sess_stim] = ...
        deal(cell(size(sessions,1),1));
    
    for isess = 1:size(sessions,1)
        
        sess = sessions(isess,1:5);
        
        load(fullfile(dirD, partid, sess, 'inv_sols_lcmv_c27_100to300ms_600to100msBl_1F2F.mat'), ...
            'lcmv_pow', 'lcmv_all', 'lcmv_pre', 'lcmv_post', ...
            'source_effect', 'source_effect_norm_sub', 'source_effect_norm_div', ...
            'source_effect_mean', 'source_effect_mean_norm_sub', 'source_effect_mean_norm_div', ...
            'W', 'eeg_post_avg', 'eeg_pre_avg');
        
        effects_lcmv{isess} = lcmv_pow.avg.pow;
        contrast_sess{isess} = source_effect;
        contrast_sess_norm_sub{isess} = source_effect_norm_sub;
        contrast_sess_norm_div{isess} = source_effect_norm_div;
        contrast_sess_mean{isess} = source_effect_mean';
        contrast_sess_mean_norm_sub{isess} = source_effect_mean_norm_sub';
        contrast_sess_mean_norm_div{isess} = source_effect_mean_norm_div';
        lcmv_pre_cell{isess} = lcmv_pre.avg.pow;
        lcmv_post_cell{isess} = lcmv_post.avg.pow;
        
        erp_pre_sess{isess} = eeg_pre_avg;
        erp_post_sess{isess} = eeg_post_avg;
    end
    
    cfg = [];
    erp_subj_pre{isubj} = ft_timelockgrandaverage(cfg, erp_pre_sess{:});
    erp_subj_post{isubj} = ft_timelockgrandaverage(cfg, erp_post_sess{:});
    
    % Average over sessions
    
    % GA effects_lcmv
    ga_lcmv{isubj} = lcmv_pre;
    tmp = nan(size(lcmv_pre.avg.pow, 1), size(sessions,1));
    for isess = 1:size(sessions,1)
        tmp(:, isess) = effects_lcmv{isess};
    end
    ga_lcmv{isubj}.avg.pow = mean(tmp,2);
    
    % GA pre lcmv
    ga_pre_lcmv{isubj} = lcmv_pre;
    tmp = nan(size(lcmv_pre.avg.pow, 1), size(sessions,1));
    for isess = 1:size(sessions,1)
        tmp(:, isess) = lcmv_pre_cell{isess};
    end
    ga_pre_lcmv{isubj}.avg.pow = mean(tmp,2);
    
    % GA post lcmv
    ga_post_lcvm{isubj} = lcmv_post;
    tmp = nan(size(lcmv_post.avg.pow, 1), size(sessions,1));
    for isess = 1:size(sessions,1)
        tmp(:, isess) = lcmv_post_cell{isess};
    end
    ga_post_lcvm{isubj}.avg.pow = mean(tmp,2);
    
    % GA contrast (over sessions) - compute for single trials, then average
    tmp = nan(size(contrast_sess{isess},2), size(sessions,1));
    for isess = 1:size(sessions,1)
        tmp(:, isess) = contrast_sess{isess}';
    end
    contrast2{isubj} = mean(tmp,2);
    
    % For each trial normalize (Post-Pre)./Pre by subtracting mean over all
    % voxels
    tmp = nan(size(contrast_sess_norm_sub{isess},2), size(sessions,1));
    for isess = 1:size(sessions,1)
        tmp(:, isess) = contrast_sess_norm_sub{isess}';
    end
    contrast2_norm_sub{isubj} = mean(tmp,2);
    
    % For each trial normalize (Post-Pre)./Pre by dividing by mean over all
    % voxels
    tmp = nan(size(contrast_sess_norm_div{isess},2), size(sessions,1));
    for isess = 1:size(sessions,1)
        tmp(:, isess) = contrast_sess_norm_div{isess}';
    end
    contrast2_norm_div{isubj} = mean(tmp,2);
    
    % GA contrast (over sessions) - average over trials, then contrast
    tmp = nan(size(contrast_sess_mean{isess},2), size(sessions,1));
    for isess = 1:size(sessions,1)
        tmp(:, isess) = contrast_sess_mean{isess}';
    end
    contrast2_mean{isubj} = mean(tmp,2);
    
    % For each session normalize (Post-Pre)./Pre by subtracting mean over all
    % voxels
    tmp = nan(size(contrast_sess_mean_norm_sub{isess},2), size(sessions,1));
    for isess = 1:size(sessions,1)
        tmp(:, isess) = contrast_sess_mean_norm_sub{isess}';
    end
    contrast2_mean_norm_sub{isubj} = mean(tmp,2);
    
    % For each session normalize (Post-Pre)./Pre by dividing by mean over all
    % voxels
    tmp = nan(size(contrast_sess_mean_norm_div{isess},2), size(sessions,1));
    for isess = 1:size(sessions,1)
        tmp(:, isess) = contrast_sess_mean_norm_div{isess}';
    end
    contrast2_mean_norm_div{isubj} = mean(tmp,2);
    
    cfg = [];
    erp_subj_pre{isubj} = ft_timelockgrandaverage(cfg, erp_pre_sess{:});
    erp_subj_post{isubj} = ft_timelockgrandaverage(cfg, erp_post_sess{:});
    
end % subj loop


%% Save / load data
mkdir(fullfile('D:\dfi_experiment_data\eeg_data\experiment\source_analysis'));
save(fullfile('D:\dfi_experiment_data\eeg_data\experiment\source_analysis', ...
    'fun_localizer_1F2F_vs_all_noise_in_coi.mat'), ...
    'leadfield', 'erp_subj_pre', 'erp_subj_post', 'contrast2*', 'lcmv_all', ...
    '-v7.3');

% Still need to fetch at the very least one lcmv_all and leadfield
load(fullfile('D:\dfi_experiment_data\eeg_data\experiment\source_analysis', ...
    'fun_localizer_1F2F_vs_all_noise_in_coi.mat'));


%% Explore source data

% Compute effect with second pipeline: 
% Source space contrast 2:
% pre_src = reshape(W * cell2mat(eeg_pre.trial), [Nvox, Ntpts, Ntrls]);
% post_src = reshape(W * cell2mat(eeg_post.trial), [Nvox, Ntpts, Ntrls]);
%  
% source_effect = mean((var(post_src,0,2)-var(pre_src,0,2))./var(pre_src,0,2),3);
% source_effect_mean = (mean(var(post_src, 0, 2),3) - ...
%                       mean(var(pre_src, 0, 2),3)) ./ ...
%                       mean(var(pre_src, 0, 2),3);

% Convert contrast cell arrays to matrices
contrast2_mat = reshape(cell2mat(contrast2), [n1, 20]);
contrast2_mat_norm_sub = reshape(cell2mat(contrast2_norm_sub), [n1, 20]);
contrast2_mat_norm_div = reshape(cell2mat(contrast2_norm_div), [n1, 20]);

contrast2_mat_mean = reshape(cell2mat(contrast2_mean), [n1, 20]);
contrast2_mat_mean_norm_sub = reshape(cell2mat(contrast2_mean_norm_sub), [n1, 20]);
contrast2_mat_mean_norm_div = reshape(cell2mat(contrast2_mean_norm_div), [n1, 20]);

% Compute grand average contrasts
GA_contrast2 = mean(reshape(cell2mat(contrast2), [n1, 20]), 2);
GA_contrast2_norm_sub = mean(reshape(cell2mat(contrast2_norm_sub), [n1, 20]), 2);
GA_contrast2_norm_div = mean(reshape(cell2mat(contrast2_norm_div), [n1, 20]), 2);

GA_contrast2_mean = mean(reshape(cell2mat(contrast2_mean), [n1, 20]), 2);
GA_contrast2_mean_norm_sub = mean(reshape(cell2mat(contrast2_mean_norm_sub), [n1, 20]), 2);
GA_contrast2_mean_norm_div = mean(reshape(cell2mat(contrast2_mean_norm_div), [n1, 20]), 2);


% Compute t-values
[tmap_ratio2, tmap_ratio2_norm_sub, tmap_ratio2_norm_div, ...
    tmap_ratio2_mean, tmap_ratio2_mean_norm_sub, tmap_ratio2_mean_norm_div] = ...
    deal(nan(size(GA_contrast2)));
disp('computing t-values over contrasts..... ')
for ivox = 1:n1
    
    [~, p, ~, stats] = ttest(contrast2_mat(ivox, :));
    tmap_ratio2(ivox) = stats.tstat;
    [~, p, ~, stats] = ttest(contrast2_mat_norm_sub(ivox, :));
    tmap_ratio2_norm_sub(ivox) = stats.tstat;
    [~, p, ~, stats] = ttest(contrast2_mat_norm_div(ivox, :));
    tmap_ratio2_norm_div(ivox) = stats.tstat;
    
    [~, p, ~, stats] = ttest(contrast2_mat_mean(ivox, :));
    tmap_ratio2_mean(ivox) = stats.tstat;
    [~, p, ~, stats] = ttest(contrast2_mat_mean_norm_sub(ivox, :));
    tmap_ratio2_mean_norm_sub(ivox) = stats.tstat;
    [~, p, ~, stats] = ttest(contrast2_mat_mean_norm_div(ivox, :));
    tmap_ratio2_mean_norm_div(ivox) = stats.tstat;
end
disp('..... done.')


%% Visualize contrasts of the form mean((post-pre)/pre)

% Plot source effect (contrast 2)
dummy = lcmv_all;
dummy.effect = nan(size(leadfield.inside));
dummy.effect(leadfield.inside) = mean(contrast2_mat, 2);

cfg = [];
cfg.parameter = 'effect';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'effect';
cfg.funcolorlim   = [min(dummy.effect) max(dummy.effect)];
ft_sourceplot(cfg, interp);
colormap('hot')
colormap('viridis')

% Plot tmap (contrast 2)
dummy.tmap = dummy.effect;
dummy.tmap(leadfield.inside) = tmap_ratio2;
cfg = [];
cfg.parameter = 'tmap';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'tmap';
cfg.funcolorlim   = [0 13];
ft_sourceplot(cfg, interp);
colormap('viridis')


% Plot source effect (contrast2 norm sub)
dummy = lcmv_all;
dummy.effect = nan(size(leadfield.inside));
dummy.effect(leadfield.inside) = mean(contrast2_mat_norm_sub, 2);
cfg = [];
cfg.parameter = 'effect';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'effect';
cfg.funcolorlim   = [min(dummy.effect) max(dummy.effect)];
ft_sourceplot(cfg, interp);
colormap('hot')
colormap('viridis')

% Plot tmap (contrast 2 norm sub)
dummy.tmap = dummy.effect;
dummy.tmap(leadfield.inside) = tmap_ratio2_norm_sub;
cfg = [];
cfg.parameter = 'tmap';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'tmap';
cfg.funcolorlim   = [0 6];
ft_sourceplot(cfg, interp);
colormap('viridis')


% Plot source effect (contrast2 norm div)
dummy = lcmv_all;
dummy.effect = nan(size(leadfield.inside));
dummy.effect(leadfield.inside) = mean(contrast2_mat_norm_div, 2);
cfg = [];
cfg.parameter = 'effect';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'effect';
cfg.funcolorlim   = [min(dummy.effect) max(dummy.effect)];
ft_sourceplot(cfg, interp);
colormap('hot')
colormap('viridis')

% Plot tmap (contrast 2 norm div)
dummy.tmap = dummy.effect;
dummy.tmap(leadfield.inside) = tmap_ratio2_norm_div;
cfg = [];
cfg.parameter = 'tmap';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'tmap';
%cfg.funcolorlim   = [0 6];
ft_sourceplot(cfg, interp);
colormap('viridis')




%% Visualize contrasts of the form (mean(post)-mean(pre))/mean(pre)

% Plot source effect (contrast 2 mean)
dummy = lcmv_all;
dummy.effect = nan(size(leadfield.inside));
dummy.effect(leadfield.inside) = mean(contrast2_mat_mean, 2);

cfg = [];
cfg.parameter = 'effect';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'effect';
cfg.funcolorlim   = [min(dummy.effect) max(dummy.effect)];
ft_sourceplot(cfg, interp);
colormap('hot')
colormap('viridis')

% Plot tmap (contrast 2 mean)
dummy.tmap = dummy.effect;
dummy.tmap(leadfield.inside) = tmap_ratio2_mean;
cfg = [];
cfg.parameter = 'tmap';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'tmap';
cfg.funcolorlim   = [-22 0];
ft_sourceplot(cfg, interp);
colormap('viridis')


% Plot source effect (contrast2 mean norm sub)
dummy = lcmv_all;
dummy.effect = nan(size(leadfield.inside));
dummy.effect(leadfield.inside) = mean(contrast2_mat_mean_norm_sub, 2);
cfg = [];
cfg.parameter = 'effect';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'effect';
cfg.funcolorlim   = [min(dummy.effect) max(dummy.effect)];
ft_sourceplot(cfg, interp);
colormap('hot')
colormap('viridis')

% Plot tmap (contrast 2 mean norm sub)
dummy.tmap = dummy.effect;
dummy.tmap(leadfield.inside) = tmap_ratio2_mean_norm_sub;
cfg = [];
cfg.parameter = 'tmap';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'tmap';
cfg.funcolorlim   = [0 6];
ft_sourceplot(cfg, interp);
colormap('viridis')


% Plot source effect (contrast2 mean norm div)
dummy = lcmv_all;
dummy.effect = nan(size(leadfield.inside));
dummy.effect(leadfield.inside) = mean(contrast2_mat_mean_norm_div, 2);
cfg = [];
cfg.parameter = 'effect';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'effect';
cfg.funcolorlim   = [min(dummy.effect) max(dummy.effect)];
ft_sourceplot(cfg, interp);
colormap('hot')
colormap('viridis')

% Plot tmap (contrast 2 mean norm div)
dummy.tmap = dummy.effect;
dummy.tmap(leadfield.inside) = tmap_ratio2_mean_norm_div;
cfg = [];
cfg.parameter = 'tmap';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'tmap';
%cfg.funcolorlim   = [-6 6];
ft_sourceplot(cfg, interp);
colormap('viridis')



% eof


