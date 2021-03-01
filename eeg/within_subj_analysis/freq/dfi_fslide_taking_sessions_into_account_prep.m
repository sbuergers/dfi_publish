% Within subject frequency sliding analysis preparation script.
% Yes-no threshold task data.
% Compute frequency sliding.
%
% Parent script(s): 
%
% Children script(s): 
%   dfi_fslide_sdtparams.m
%   dfi_fslide_taking_sessions_into_account_sdtparams.m
%
% Sibling script(s):
%   dfi_fslide_yesno_taking_sessions_into_account_prep.m
%
%
% DETAILS
%
% Compute frequency sliding for -1.2 to 0.7s time window relative
% to first stimulus onset for each subject for channels
% O1, Oz, O2, PO8, PO4, POz.
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
an_fold = 'inst_freq_ynt';

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
catch
    warning('Cannot find fieldtrip folder')
end

subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N        = 20;
task     = 'yn_threshold';
eegfile  = 'data_preproc2500_to4500.mat';

% Initialize variables for iAF estimation
Fs         = 256;      % Downsample to this frequency
w          = [6 14];   % alpha search window


for isubj = 1:N
        
    subj_dir = fullfile(data_dir, subjvect{isubj}, task);
    
    % get participant, session and session/run information from user
    partid    = subjvect{isubj};
    rec_type  = 'eyes_open';
    resp_lock = 'no';
    stim_lock = 'yes';
    
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
                
                % select trial window (otherwise we run out of memory)
                cfg          = [];
                cfg.toilim   = [-1.2, 0.7];  % 0.4 s post-stim are later set to zero
                eeg_cell{cnt} = ft_redefinetrial(cfg, data_preproc);
                clear data_preproc
                
                sessions = [sessions, isess];
                runs     = [runs, irun];
                if strcmp(rec_type, 'eyes_open')
                    beh_cell{cnt} = bdata; clear bdata
                    load(fullfile(temp_dir, 'artifact_info.mat'));
                    art_cell{cnt} = artifact_samples; clear artifact_samples
                    if exist('eye_behmatch', 'var')
                        eye_cell{cnt} = eye_behmatch; clear eye_behmatch
                        Iev_cell{cnt} = data_eye; clear data_eye
                    end
                end
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
    
    
    %% 1.) --- PREPROCESSING ---
    
    for ifile = 1:numel(eeg_cell)
        
        fprintf('\n\nPreprocessing stimulus locked data...\n\n');
        
        % get behavioral and artifact data of this recording
        bdata            = beh_cell{ifile};
        artifact_samples = art_cell{ifile};

        % delete trials that contain artifacts
        cfg                           = [];
        cfg.artfctdef.muscle.artifact = artifact_samples;
        cfg.artfctdef.reject          = 'complete';
        
        stim_clean = ft_rejectartifact(cfg, eeg_cell{ifile});
        
        % adjust behavioral data accordingly
        oldsampinfo = eeg_cell{ifile}.sampleinfo(:,1);
        newsampinfo = stim_clean.sampleinfo(:,1);
        keepTrials  = find( ismember(oldsampinfo, newsampinfo) );
        
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
        
        % Downsample and demean
        cfg            = [];
        cfg.resamplefs = Fs;
        cfg.detrend    = 'no';
        fprintf('\n\nDownsampling to %i Hz...\n\n', cfg.resamplefs);
        stim_downsamp  = ft_resampledata(cfg, stim_clean);
        
        % detrend
        cfg          = [];
        cfg.detrend  = 'yes';
        cfg.overlap  = 0;
        stim_preproc = ft_preprocessing(cfg, stim_downsamp);
        
        % Rereference
        fprintf('\n\nRereferencing...\n\n');
        cfg            = [];
        cfg.channel    = {'all'; '-blink'; '-pupil'; '-gazeX'; '-gazeY'};
        cfg.reref      = 'yes';
        cfg.refchannel = 'all';
        stim_reref     = ft_preprocessing(cfg, stim_preproc);
        
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
        stim_reref   = ft_preprocessing(cfg, stim_reref);
        
        eeg_cell_stim{ifile} = stim_reref;
        beh_cell_stim{ifile} = bdata_stim;
        
        clear stim_downsamp stim_reref stim_preproc stim_seg
        
    end % file loop
    
    clear eeg_cell eye_cell
    
    
    %
    %% 2.) --- SESSION SELECTION ---
    %
    ncomb = length(unique(sessions))+1;
    sessionIDs = unique(sessions);
    nchan = nan(1, length(unique(sessions)));
    
    for isess_combi = ncomb
        
        if isess_combi == ncomb % all sessions combined
            
            % folders
            sess_fold = 'all_yn_sessions';
            
            % append data
            if numel(eeg_cell_stim) > 1
                cfg      = [];
                eeg_stim = ft_appenddata(cfg, eeg_cell_stim{:});
            else
                eeg_stim = eeg_cell_stim{1};
            end
            
            beh_stim = [];
            for i = 1:numel(beh_cell_stim)
                beh_stim = [beh_stim; beh_cell_stim{i}];
            end
            
        else % single sessions
            
            % folders
            sess_fold = sprintf('session_%i', sessionIDs(isess_combi));
            
            % append data
            if sum(sessions==sessionIDs(isess_combi)) > 1
                cfg      = [];
                eeg_stim = ft_appenddata(cfg, eeg_cell_stim{sessions==sessionIDs(isess_combi)});
            else
                eeg_stim = eeg_cell_stim{sessions==sessionIDs(isess_combi)};
            end
            
            beh_stim = [];
            ids_of_this_session = find(sessions==sessionIDs(isess_combi));
            for i = ids_of_this_session
                beh_stim = [beh_stim; beh_cell_stim{i}];
            end
            
        end
        
        % do eeg and bdata still match?
        disp('Checking if behavioral and EEG data still match...')
        behid   = beh_stim.trlid;
        eegid   = eeg_stim.trialinfo;
        id_diff = eegid - behid;
        if any(mod(id_diff, 10)) ~= 0
            [behid, eegid]
            error('eeg and behavioral data do not match!');
        else
            disp('Phew, data still match!')
        end
        
        % Get trial IDs
        [id_stim, idm_stim] = dfi_get_trial_indices( beh_stim, 'yes' );
   


        %% 4.) --- FREQUENCY SLIDING ---
   
        s = dfi_eeg_defaults;
        s.chan      = {'O*', 'POz', 'PO4', 'PO8'};
        s.psd       = 0;
        s.tfr       = 0;
        s.plv       = 0;
        s.phase     = 0;
        s.fslide    = 1;
        s.filt      = 'no';    % no initial bp filter, we use a tukeyfilter anyway
        s.alpha     = w;
        s.foi       = s.alpha;
        s.tukeyfilt = 'firls'; % Cohen uses 'firls' with 15% transition,
        % fieldtrip uses 'firws' with half-attenuation at
        % foi(1) and foi(2)
        s.fsplot    = 0;       % Plot processing steps?
        
		% select time window of interest + padding
        cfg         = [];
        cfg.toilim  = [-1.2 0];
        eeg_pre  = ft_redefinetrial(cfg, eeg_stim);
        
        % Zero pad post stim period
        s.padwin    = [-1.2 0.7];
        eeg_zeropad = eeg_pre;
        for itrl = 1:numel(eeg_pre.trial)
            eeg_zeropad.time{itrl}  = [eeg_pre.time{itrl}, -eeg_pre.time{itrl}(end-1:-1:ceil(length(eeg_pre.time{itrl})/2))];
            eeg_zeropad.trial{itrl} = [eeg_pre.trial{itrl}(:,:), zeros(64,ceil(length(eeg_pre.time{itrl})/2))];
        end
		assert(size(eeg_zeropad.trial{1}, 2) == size(eeg_zeropad.time{1}, 2))

        % Compute frequency sliding for all trials
        [~, eeg_fslide]  = dfi_eeg_analysis(s, eeg_zeropad);     

        % Save data
        mkdir(fullfile(subj_dir, an_fold, sess_fold));
        save(fullfile(subj_dir, an_fold, sess_fold, 'fslide_data_trls_all'), ...
            'eeg_fslide', 'beh_stim', '-v7.3');

        % Clean up a bit
        close all
        clear eeg_fslide eeg_zeropad eeg_pre
        
   
    end % session combinations
    
    clear eeg_cell_stim beh_cell_stim
    
end % loop over participants


% enable warnings again
warning('on','all')


% // eof

















