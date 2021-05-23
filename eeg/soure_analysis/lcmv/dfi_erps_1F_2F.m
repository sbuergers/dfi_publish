% Compute scalp and source ERPs for 1F and 2F trials.
%
% Parent script(s): 
%   dfi_source_proj_erp.m
%
% Children script(s): 
%   None
%
% Sibling script(s):
%   None
%
% DETAILS
% Compute scalp and source ERPs for 1F and 2F trials. This is done
% separately from the main source analysis scripts, because there we used
% a discontinuous chunk of data (-0.6 to -0.1s pre-stim, 0.1 to 0.3s post-
% stim). In addition we pooled over 1F and 2F trials. It might instead be
% nice to see 1F and 2F ERPs continously from -0.6 to 0.3 s separately for
% 1F and 2F trials.
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
            
            % Get filter weights of source analysis computed in
            % dfi_source_proj_erp.m.
            load(fullfile(src_dir, partid, sess, ...
                'inv_sols_lcmv_c27_100to300ms_600to100msBl_1F2F.mat'), ...
                'W', 'leadfield');

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
            
            % Get trial IDs
            if length(unique(beh_stim.resp)) == 2
                [id_stim, idm_stim] = dfi_get_trial_indices( beh_stim, 'no' );
            else
                [id_stim, idm_stim] = dfi_get_trial_indices( beh_stim, 'yes' );
            end
            
            % EEG data (1F and 2F trials)
            cfg = [];
            cfg.trials = id_stim.v1;
            eeg_1f = ft_selectdata(cfg, eeg_stim);
            cfg.trials = id_stim.v2;
            eeg_2f = ft_selectdata(cfg, eeg_stim);
            
            % Baseline correct
            cfg = [];
            cfg.demean = 'yes';
            cfg.baselinewindow = [-0.6 -0.1];
            eeg_1f = ft_preprocessing(cfg, eeg_1f);
            eeg_2f = ft_preprocessing(cfg, eeg_2f);
            
            % pre-stim and post-stim (noise and signal)
            cfg = [];
            cfg.toilim = [-0.6, +0.3];
            eeg_1f = ft_redefinetrial(cfg, eeg_1f);
            eeg_2f = ft_redefinetrial(cfg, eeg_2f);
            
            % compute the covariance matrix for AVG
            cfg = [];
            eeg_1f_trls_avg = ft_timelockanalysis(cfg, eeg_1f);
            eeg_2f_trls_avg = ft_timelockanalysis(cfg, eeg_2f);
            
            % Also compute contrasts by hand (more control)
            Nvox = sum(leadfield.inside);
            Ntpts = size(eeg_1f.time{1}, 2);
            Ntrls = length(eeg_1f.trialinfo);
            eeg_1f_src = reshape(W * cell2mat(eeg_1f.trial), [Nvox, Ntpts, Ntrls]);
            
            Ntrls = length(eeg_2f.trialinfo);
            eeg_2f_src = reshape(W * cell2mat(eeg_2f.trial), [Nvox, Ntpts, Ntrls]);
            
            % Plot source and scalp ERPs
            load(fullfile(src_dir, 'ROI_1f2f_vs_all_noise_in_coi.mat'), ...
                'lf_vox_for_centroid', 'centroids');
            voxoi = nan(size(lf_vox_for_centroid{1}, 1), 1);
            for ivox = 1:size(lf_vox_for_centroid{1}, 1)
                voxoi(ivox) = find(ismember(leadfield.pos, ...
                    lf_vox_for_centroid{1}(ivox, :),'rows'));
            end
            
            figure;
            time = eeg_1f.time{1};
            choi = find(ismember(eeg_1f.label, 'O2'));
            eeg_1f_src_avg = squeeze(nanmean(nanmean(eeg_1f_src(voxoi,:,:)), 3));
            eeg_2f_src_avg = squeeze(nanmean(nanmean(eeg_2f_src(voxoi,:,:)), 3));
            subplot(221)
            plot(time, eeg_1f_trls_avg.avg(choi, :));
            ylabel('Sensor level')
            title('1 flash trials')
            subplot(222)
            plot(time, eeg_2f_trls_avg.avg(choi, :));
            title('2 flash trials')
            subplot(223)
            plot(time, eeg_1f_src_avg);
            ylabel('Source level')
            subplot(224)
            plot(time, eeg_2f_src_avg);
            close all
            
            % save filters and contrasts
            mkdir(fullfile(src_dir, partid, sess))
            save(fullfile(src_dir, partid, sess, ...
                'erps_min600to300ms_1F2F.mat'), ...
                'eeg_1f_src', 'eeg_2f_src', 'eeg_1f', 'eeg_2f', ...
                'eeg_1f_trls_avg', 'eeg_2f_trls_avg', 'voxoi', ...
                '-v7.3');
            
            clear eeg_trls source eeg_stim
            
        end  % session combinations
        
        clear eeg_cell beh_cell
        
    end  % loop over participants
    
end % loop over tasks

% enable warnings again
warning('on','all')


% eof
