% Alpha peak frequency estimation for pre-stimulus, task related data
% in source space (both yes-no and yes-no threshold tasks).
% Use alpha peak frequency fitting method by Corcoran et al. (2018)
%
% Dependencies:
% 	restingIAF (Corcoran et al. toolbox), with modified scripts
% 	meanIAF_sb.m
% 	restingIAF_sb.m
%
% Parent script(s): 
%   dfi_source_proj_erp.m
%
% Children script(s): 
%   ?
%
% Sibling script(s):
%   dfi_iAF_fits_eyes_open.m
%   dfi_iAF_fits_eyes_open_ynt.m
%	dfi_iAF_fits_eyes_closed_w_0padding.m
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
tasks    = {'yesno', 'yn_threshold'};

% Get voxels of interest (ROI)
load(fullfile(src_dir, 'ROI_1f2f_vs_all_noise_in_coi.mat'), ...
    'lf_vox_for_centroid', 'centroids');
roi = lf_vox_for_centroid{1};

% Get a leadfield structure (for .pos field to find ROI)
load(fullfile(src_dir, '701', 'sess2', ...
	'inv_sols_lcmv_c27_100to300ms_600to100msBl_1F2F.mat'), 'leadfield');
pos_inside = leadfield.pos(leadfield.inside,:);
roi_idx = ismember(pos_inside, roi, 'rows');


% Initialize variables for iAF estimation
Fs         = 256; % Downsample to this frequency
[pSpec]    = deal(struct('chans', [], 'sums', [])); % Contains fits for each 
                                                    % subejct and session.
                                                    % The last session
                                                    % field is the average
                                                    % over all sessions
muPaf      = nan(N, 1); % Here I average over peak estimates for all sessions
muCog      = nan(N, 1);
muPow      = nan(N, 1);
muPow_osc     = nan(N, 1);
muPow_pk      = nan(N, 1);
muPow_pk_osc  = nan(N, 1);
sessionmat = nan(N, 7);

% Peak fitting parameter settings
cmin   = 3;        % min # of channels required for cross-channel average
fRange = [3 30];   % frequency range going into fitting process
w      = [6 14];   % alpha peak search window
Fw     = 11;       % Window of polinomial filter 
k      = 5;        % Order of polinomial filter (SGT)
mpow   = 1;        % error bound used to determine threshold differentiating 
                   % substantive peaks from background spectral noise 
mdiff  = 0.2;      % minimum difference between peaks
tlen   = ceil(0.5*Fs); % window length (corresponds to 0.5 s)
tover  = 0;        % overlap between windows for welch
nfft   = 8*Fs;     % length to zero pad the data to
norm   = true;     % normalize spectra by average power


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
            
            % pre-stim data
            cfg = [];
            cfg.trials = 'all';
            eeg_pre = ft_selectdata(cfg, eeg_stim);
            
            % Baseline correct
            cfg = [];
            cfg.demean = 'yes';
            cfg.baselinewindow = [-0.6 -0.1];
            eeg_pre = ft_preprocessing(cfg, eeg_pre);
            
            % pre-stim and post-stim (noise and signal)
            cfg = [];
            cfg.toilim = [-0.6, -0.1];
            eeg_pre = ft_redefinetrial(cfg, eeg_pre);
                        
            % Load source filters
            load(fullfile(src_dir, partid, sess, ...
                'inv_sols_lcmv_c27_100to300ms_600to100msBl_1F2F.mat'), 'W');
                        
            % Project data into source space
            Ntpts = size(eeg_pre.time{1}, 2);
            Ntrls = length(eeg_pre.trialinfo);
            Nvchan = sum(roi_idx);
            
            src_oi = nan(Nvchan, Ntpts, Ntrls);
            for itrl=1:Ntrls
                src = W * eeg_pre.trial{itrl};
                src_oi(:, :, itrl) = src(roi_idx,:,:);
            end

            % Create fake ft structure
            vchan_labels = cell(Nvchan, 1);
            for ivchan = 1:Nvchan
                vchan_labels{ivchan} = sprintf('src%i', ivchan);
            end
            eeg_src = eeg_stim;
            eeg_src.label = vchan_labels;
            for itrl = 1:size(src_oi,3)
                eeg_src.trial{itrl} = src_oi(:,:,itrl);
            end            
            
            %
            %% 3.) --- ALPHA PEAK FREQUENCY ---
            %
            jx = isess;

            fprintf('\n\nCorcoran iAF estimation in source space:')
            fprintf('Subject: %i, session %i\n\n', isubj, sessions(jx));

            % treat trials as artifically continous segment (pwelch is going to cut it
            % up again exactly in the same way, and here we use no overlap!)
            ps_trls = eeg_src;
            lentrl = size(ps_trls.trial{1},2);
            fk_cont_data = nan(numel(ps_trls.label), numel(ps_trls.trial)*lentrl);
            for i = 1:numel(ps_trls.trial)
                fk_cont_data(:, (1+(i-1)*lentrl):(lentrl*i) ) = single(ps_trls.trial{i});
            end

            nchan(jx) = length(ps_trls.label);

            [pSpec(isubj, jx).sums, pSpec(isubj, jx).chans, f] = restingIAF_sb(...
                fk_cont_data, nchan(jx), cmin, fRange, Fs, w, Fw, ...
                k, mpow, mdiff, tlen, tover, nfft, norm ...
            );

        end  % session combinations
        
        % weighted average of mean IAF estimates across (j-th) recordings
        if length(unique(sessions)) > 1
            [muPaf(isubj,:), muCog(isubj,:), muPow(isubj,:), ...
             muPow_osc(isubj,:), muPow_pk(isubj,:), muPow_pk_osc(isubj,:)] = ...
                meanIAF_sb([pSpec(isubj, 1:length(unique(sessions))).sums], ...
                    nchan(1:length(unique(sessions))) ...
                );
        end
                
    end  % loop over participants
    
    figure
    plot(muPaf, 'b'); hold on
    title(sprintf('Individual alpha peak frequencies (%s)', task))
    xlabel('Subject ID'); set(gca, 'xtick', 1:20);
    ylabel('iAF (Hz)'); grid on; xlim([1 20])
    
    
    %% 4.) --- SAVE DATA ---
    
    if ~exist(fullfile(fig_dir, an_fold), 'dir')
        mkdir(fullfile(fig_dir, an_fold));
    end
    
    fname = 'eyes_open_pkinfo';
    mkdir(fullfile(src_dir, task));
    save(fullfile(src_dir, task, fname), ...
        'muPaf*', 'muPow*', 'pSpec', 'cmin', 'channel_vect', ...
        'fRange', 'w', 'Fs', 'Fw', 'k');
    
end % loop over tasks


% enable warnings again
warning('on','all')


% eof

