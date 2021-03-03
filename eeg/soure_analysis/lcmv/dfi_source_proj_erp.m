% Compute spatial filters for source projection based on ERPs
%
% Parent script(s): 
%   None
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
% Compute filters for evenly spaced grid in source space
% using linearly constrained minimum variance beamforming
% following Van Veen et al. (1997). 
% 
% We use the -0.6 to -0.1s pre-stimulus period and 0.1 to 0.3s
% post-stimulus period to compute the covariance matrix necessary
% for estimating the filters. For the pre-stimulus data we include
% all conditions (as all conditions are used for the analyses in
% source space later on), but for post-stim we only focus on
% 1-flash and 2-flash conditions (i.e. without beeps), to avoid
% caveats related to correlated sources. 
%
% To assess filter quality, we compute various contrasts of the
% form (VAR(post)-VAR(pre) ./ VAR(pre) in source space. Specifically, 
% we either compute the contrast for single trials and then average,
% or average over trials and then compute the contrast. Moreover,
% we try two normalizations: 1. Subtract average contrast value over
% all voxels, and 2. Divide by average contrast value over all
% voxels. 
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


restoredefaultpath; clc; close all; clear all

% disable warnings for speed sake:
warning('off','all')

% Initialization
try
    addpath(genpath('dfi'))
    dfi_startup
catch
    warning('Cannot find dfi folder')
end
ft_path = fullfile('fieldtrip20200906', 'fieldtrip-master'); 
addpath(ft_path);
ft_defaults

data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
src_dir = fullfile(data_dir, 'source_analysis');
eegfile  = 'data_preproc2500_to4500.mat';

subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N        = length(subjvect);
Fs       = 64;  % Downsample to this frequency
tasks    = {'yesno', 'yn_threshold'};

ft_headmodel_folder = fullfile(ft_path, 'template', 'headmodel');
mri = ft_read_mri(fullfile(ft_headmodel_folder, 'standard_mri.mat'));
seg = ft_read_mri(fullfile(ft_headmodel_folder, 'standard_seg.mat'));
headmodel = load(fullfile(ft_headmodel_folder, 'standard_bem.mat'));
headmodel = headmodel.vol;

cfg             = [];
cfg.resolution  = 8;  % in mm
cfg.headmodel   = headmodel;
cfg.inwardshift = 1;  % shifts dipoles away from surfaces
sourcemodel     = ft_prepare_sourcemodel(cfg);


% Loop over tasks, subjects and sessions
for itask = 1:2
    task  = tasks{itask};
    
    for isubj = 1:N
        
        %
        %% 1.) ---- Preprocessing ----
        %
        fprintf('\n\n\n\n\n\nTASK %s, PARTICIPANT %s\n\n\n\n\n\n', task, subjvect{isubj})
        
        subj_dir = fullfile(data_dir, subjvect{isubj}, task);
        
        % get participant, session and session/run information from user
        partid    = subjvect{isubj};
        rec_type  = 'eyes_open';
        
        % get all data files and append
        cnt = 1;
        sessions = [];
        runs     = [];
        for isess = 1:10
            for irun = 1:10
                temp_dir = fullfile(subj_dir, sprintf('sess%d', isess),  ...
                    sprintf('run%d', irun), rec_type);
                if exist(fullfile(temp_dir, eegfile), 'file')
                    fprintf('Loading data from \n%s\n\n', temp_dir)
                    load(fullfile(temp_dir, eegfile));
                    sessions = [sessions, isess];
                    runs     = [runs, irun];
                    if strcmp(rec_type, 'eyes_open')
                        load(fullfile(temp_dir, 'artifact_info.mat'));
                    end
                    
                    % select trial window
                    cfg          = [];
                    cfg.toilim   = [-0.6, 0.6];
                    stim_seg     = ft_redefinetrial(cfg, data_preproc);
                    clear data_preproc
                    
                    % delete trials that contain artifacts
                    cfg                           = [];
                    cfg.artfctdef.muscle.artifact = artifact_samples;
                    cfg.artfctdef.reject          = 'complete';
                    
                    stim_clean = ft_rejectartifact(cfg, stim_seg);
                    
                    % adjust behavioral data accordingly
                    oldsampinfo = stim_seg.sampleinfo(:,1);
                    newsampinfo = stim_clean.sampleinfo(:,1);
                    keepTrials  = find( ismember(oldsampinfo, newsampinfo) );
                    clear stim_seg
                    
                    bdata_stim = bdata(keepTrials, :);
                    
                    % do eeg and bdata still match?
                    disp('Checking if behavioral and EEG data still match...')
                    behid   = bdata_stim.trlid;
                    eegid   = stim_clean.trialinfo;
                    id_diff = eegid - behid;
                    if any(mod(id_diff, 10)) ~= 0
                        [behid, eegid]
                        error('eeg and behavioral data do not match!');
                    else
                        disp('Phew, data still match!')
                    end
                    
                    % Filter
                    cfg              = [];
                    cfg.lpfilter     = 'yes';
                    cfg.lpfreq       = 28;
                    cfg.lpfiltord    = 10;
                    cfg.lpfiltdir    = 'twopass';
                    data_lowpass     = ft_preprocessing(cfg, stim_clean);
                    clear stim_clean
                    
                    % Downsample
                    cfg            = [];
                    cfg.resamplefs = Fs;
                    cfg.detrend    = 'no';
                    fprintf('\n\nDownsampling to %i Hz...\n\n', cfg.resamplefs);
                    stim_downsamp  = ft_resampledata(cfg, data_lowpass);
                    clear data_lowpass
                    
                    % Rereference
                    fprintf('\n\nRereferencing...\n\n');
                    cfg            = [];
                    cfg.channel    = {'all'; '-blink'; '-pupil'; '-gazeX'; '-gazeY'};
                    cfg.reref      = 'yes';
                    cfg.refchannel = 'all';
                    stim_reref     = ft_preprocessing(cfg, stim_downsamp);
                    clear stim_downsamp
                    
                    % delete trials that contain NaNs for some reason
                    nanTrialIds = [];
                    for i=1:size(stim_reref.trialinfo,1)
                        if sum(any(isnan(stim_reref.trial{i}))) > 0
                            nanTrialIds = [nanTrialIds; i];
                        end
                    end
                    
                    sprintf('\n\n The following trials contain NaNs for some reason... delete them!')
                    nanTrialIds
                    bdata_stim(nanTrialIds,:) = [];
                    
                    cfg          = [];
                    cfg.trials   = setdiff(1:size(stim_reref.trialinfo,1), nanTrialIds);
                    data_preproc = ft_preprocessing(cfg, stim_reref);
                    
                    % combine runs and sessions in cell array
                    if strcmp(rec_type, 'eyes_open')
                        beh_cell{cnt} = bdata_stim; clear bdata
                    end
                    eeg_cell{cnt} = data_preproc;  clear data_preproc
                    cnt = cnt + 1;
                end
            end
        end
        
        if ~exist('eeg_cell', 'var')
            error('Ups. Could not find EEG data ... Change path so that you can see dfi_experiment_data.')
        end
        
        % session name
        session = 'all_sessions';
        supertitle = ['Participant ', partid, ', all sessions', ', all runs'];
        
        %
        %% 2.) --- Source projection per session ---
        %
        ncomb = length(unique(sessions));
        sessionIDs = unique(sessions);
        nchan = nan(1, length(unique(sessions)));
        
        for isess = 1:ncomb
            
            sess = sprintf('sess%i', sessionIDs(isess));

            % combine EEG and behavioural data for all runs of this sessions
            if sum(sessions==sessionIDs(isess)) > 1
                cfg      = [];
                eeg_stim = ft_appenddata(cfg, eeg_cell{sessions==sessionIDs(isess)});
            else
                eeg_stim = eeg_cell{sessions==sessionIDs(isess)};
            end
            
            beh_stim = [];
            ids_of_this_session = find(sessions==sessionIDs(isess));
            for i = ids_of_this_session
                beh_stim = [beh_stim; beh_cell{i}];
            end
            
            % do eeg and bdata still match?
            disp('Checking if behavioral and EEG data still match...')
            behid   = beh_stim.trlid;
            eegid   = eeg_stim.trialinfo;
            id_diff = eegid - behid;
            if any(mod(id_diff, 10) ~= 0)
                [behid, eegid]
                error('eeg and behavioral data do not match!');
            else
                disp('Phew, data still match!')
            end
            
            % Load electrode positions aligned to head surface in MNI space
            load(fullfile(src_dir, partid, sess, 'elec_aligned.mat'))
            
            % project electrodes onto scalp surface
            cfg = [];
            cfg.headshape = headmodel.bnd(1);
            cfg.method = 'project';
            [elec_surf] = ft_electroderealign(cfg, elec_aligned);
            
            % compute leadfield
            cfg                  = [];
            cfg.headmodel        = headmodel;  % volume conduction headmodel
            cfg.sourcemodel      = sourcemodel;  % normalized grid positions
            cfg.elec             = elec_surf;
            leadfield            = ft_prepare_leadfield(cfg);
            
            % Get trial IDs
            if length(unique(beh_stim.resp)) == 2
                [id_stim, idm_stim] = dfi_get_trial_indices( beh_stim, 'no' );
            else
                [id_stim, idm_stim] = dfi_get_trial_indices( beh_stim, 'yes' );
            end
            
            % pre-stim data
            cfg = [];
            cfg.trials = 'all';
            eeg_pre = ft_selectdata(cfg, eeg_stim);
            
            % Post-stim data (1F and 2F trials)
            cfg = [];
            cfg.trials = [id_stim.v1; id_stim.v2];
            eeg = ft_selectdata(cfg, eeg_stim);
            
            % Baseline correct
            cfg = [];
            cfg.demean = 'yes';
            cfg.baselinewindow = [-0.6 -0.1];
            eeg_post = ft_preprocessing(cfg, eeg);
            eeg_pre = ft_preprocessing(cfg, eeg_pre);
            
            % pre-stim and post-stim (noise and signal)
            cfg = [];
            cfg.toilim = [-0.6, -0.1];
            eeg_pre = ft_redefinetrial(cfg, eeg_pre);
            
            cfg = [];
            cfg.toilim = [+0.1, +0.3];
            eeg_post = ft_redefinetrial(cfg, eeg_post);
            
            % combine data in one structure
            cfg = [];
            eeg_comb = ft_appenddata(cfg, eeg_pre, eeg_post);
            
            % compute the covariance matrix for AVG
            cfg = [];
            cfg.covariance = 'yes';
            cfg.removemean = 'yes';  % default = 'yes'
            eeg_trls_avg = ft_timelockanalysis(cfg, eeg_comb);
            eeg_pre_avg = ft_timelockanalysis(cfg, eeg_pre);
            eeg_post_avg = ft_timelockanalysis(cfg, eeg_post);
            
            % Fix covariance computation for eeg_trls_avg (only necessary when
            % there is a time gap between pre and post).
            trl_covs = zeros(64, 64);
            nsamp = 0;
            for itrl = 1:size(eeg_comb.trial,2)
                dat = eeg_comb.trial{itrl} - repmat(mean(eeg_comb.trial{itrl}), [64, 1]);
                trl_covs = trl_covs + dat * dat';
                nsamp = nsamp + size(eeg_comb.trial{itrl}, 2) - 1;
            end
            sum(sum(abs(trl_covs / nsamp - eeg_trls_avg.cov)))
            
            % Use my cov
            eeg_trls_avg.cov = trl_covs / nsamp;
            
            % LCMV
            cfg                 = [];
            cfg.method          = 'lcmv';
            cfg.sourcemodel     = leadfield;
            cfg.elec            = elec_surf;
            cfg.senstype        = 'EEG';
            cfg.keepleadfield   = 'yes';
            cfg.lcmv.keepfilter = 'yes';
            cfg.lcmv.fixedori   = 'yes';
            cfg.lcmv.lambda     = '5%';
            lcmv_all = ft_sourceanalysis(cfg, eeg_trls_avg);
            
            cfg.sourcemodel.filter = lcmv_all.avg.filter;
            lcmv_pre = ft_sourceanalysis(cfg, eeg_pre_avg);
            lcmv_post = ft_sourceanalysis(cfg, eeg_post_avg);
            
            lcmv_pow = lcmv_pre;
            lcmv_pow.avg.pow = (lcmv_post.avg.pow-lcmv_pre.avg.pow)./lcmv_pre.avg.pow;
            
            %                 % plot
            %                 dummy = lcmv_all;
            %                 dummy.effect = lcmv_pow.avg.pow;
            %
            %                 cfg = [];
            %                 cfg.parameter = 'effect';
            %                 [interp] = ft_sourceinterpolate(cfg, dummy, mri);
            %
            %                 cfg = [];
            %                 cfg.method        = 'slice';
            %                 cfg.funparameter  = 'effect';
            %                 cfg.funcolorlim   = [-max(dummy.effect) max(dummy.effect)];
            %                 ft_sourceplot(cfg, interp);
            %                 colormap('viridis')
            
            
            % Also compute contrasts by hand (more control)
            Nvox = sum(leadfield.inside);
            Ntpts = size(eeg_pre.time{1}, 2);
            Ntrls = length(eeg_pre.trialinfo);
            
            W = cell2mat(lcmv_all.avg.filter(leadfield.inside));
            
            pre_src = reshape(W * cell2mat(eeg_pre.trial), [Nvox, Ntpts, Ntrls]);
            
            Ntpts = size(eeg_post.time{1}, 2);
            Ntrls = length(eeg_post.trialinfo);
            post_src = reshape(W * cell2mat(eeg_post.trial), [Nvox, Ntpts, Ntrls]);
            
            % get source effect estimate
            Nrand = 100;
            [src_eff, src_eff_norm_sub, src_eff_norm_div] = deal(nan(Nrand, sum(leadfield.inside)));
            for i = 1:Nrand
                fprintf('randomization %i\n', i)
                trlsel = randsample(length(eeg_pre.trialinfo), Ntrls);
                src_eff(i,:) = mean( (var(post_src, 0, 2) - ...
                    var(pre_src(:,:,trlsel), 0, 2)) ./ ...
                    var(pre_src(:,:,trlsel), 0, 2), 3);
                
                % center each trial's contrast before averaging over trials
                X = (var(post_src, 0, 2) - ...
                    var(pre_src(:,:,trlsel), 0, 2)) ./ ...
                    var(pre_src(:,:,trlsel), 0, 2);
                X_vox_avg = squeeze(mean(X, 1));
                src_eff_norm_sub(i,:) = mean( squeeze(X) - ...
                    repmat(X_vox_avg', [size(X, 1), 1]) ...
                    , 2);
                src_eff_norm_div(i,:) = mean( squeeze(X) ./ ...
                    repmat(X_vox_avg', [size(X, 1), 1]) ...
                    , 2);
            end
            source_effect = mean(src_eff);
            source_effect_norm_sub = mean(src_eff_norm_sub);
            source_effect_norm_div = mean(src_eff_norm_div);
            
            
            %                                         % plot
            %                                         dummy = lcmv_all;
            %                                         dummy.effect = nan(size(leadfield.inside));
            %                                         dummy.effect(leadfield.inside) = source_effect;
            %
            %                                         cfg = [];
            %                                         cfg.parameter = 'effect';
            %                                         [interp] = ft_sourceinterpolate(cfg, dummy, mri);
            %
            %                                         cfg = [];
            %                                         cfg.method        = 'slice';
            %                                         cfg.funparameter  = 'effect';
            %                                         cfg.funcolorlim   = [-max(dummy.effect) max(dummy.effect)];
            %                                         ft_sourceplot(cfg, interp);
            %                                         colormap('viridis')
            
            % First average over variance for pre and post-stim, then contrast
            source_effect_mean = (mean(var(post_src, 0, 2),3) - ...
                mean(var(pre_src, 0, 2),3)) ./ ...
                mean(var(pre_src, 0, 2),3);
            
            % normalize mean source effect subtracting voxel average
            source_effect_mean_norm_sub = source_effect_mean - mean(source_effect_mean);
            
            % normalize mean source effect dividing by voxel average
            source_effect_mean_norm_div = source_effect_mean ./ mean(source_effect_mean);
            
            % Compute contrast at the scalp for comparison
            Nchan = length(eeg_post.label);
            Ntpts = length(eeg_pre.time{1});
            Ntrls = length(eeg_pre.time);
            pre_scalp = reshape(cell2mat(eeg_pre.trial), [Nchan, Ntpts, Ntrls]);
            
            Ntpts = length(eeg_post.time{1});
            Ntrls = length(eeg_post.time);
            post_scalp = reshape(cell2mat(eeg_post.trial), [Nchan, Ntpts, Ntrls]);
            
            scalp_effect_mean = (mean(var(post_scalp, 0, 2),3) - ...
                                 mean(var(pre_scalp, 0, 2),3)) ./ ...
                                 mean(var(pre_scalp, 0, 2),3);

            % save filters and contrasts
            mkdir(fullfile(src_dir, partid, sess))
            save(fullfile(src_dir, partid, sess, 'inv_sols_lcmv_c27_100to300ms_600to100msBl_1F2F.mat'), ...
                'source_effect*', 'W', 'lcmv*', 'scalp_effect_mean', ...
                'eeg_post_avg*', 'eeg_pre_avg*', 'leadfield', ...
                '-v7.3');
            
            clear eeg_trls source eeg_stim
            
        end  % session combinations
        
        clear eeg_cell beh_cell
        
    end  % loop over participants
    
end % loop over tasks

% enable warnings again
warning('on','all')


% eof
