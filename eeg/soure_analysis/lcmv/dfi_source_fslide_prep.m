function dfi_source_fslide_prep(task, roi, ft_path, data_dir)
% dfi_source_fslide(task, roi, ft_path, data_dir)
%
% DESCRIPTION
% Computes instantaneous frequency from -0.6 to -0.1s pre-stimulus
% on source projected data based on filters from LCMV ERP source localization.
%
% INPUT
% task:     Either 'yesno' or 'yn_threshold'
% roi:      N x 3 matrix of leadfield voxel coordinates for which to project data
% ft_path:  Fieldtrip main path (default = 'fieldtrip20200906\fieldtrip-master')
% data_dir: Main data directory (default = 'dfi_experiment_data\eeg_data\experiment')
%
%
% Parent script(s): 
%   dfi_source_determine_func_roi.m
%
% Children script(s): 
%   ?.m
%
% Sibling script(s):
%   None
%
%
% DETAILS
%
% Compute frequency sliding for -1.2 to 0.7s time window relative
% to first stimulus onset for each subject and session in source space
% based on LCMV filters from dfi_source_proj_erp.m for voxels specified by roi. 
%
% NOTE
% Even though this is a function, several paths are still hard coded.
% They need to be modified manually depending on where data and scripts are 
% located.
% 
% ---
% Steffen Buergers, sbuergers@gmail.com,
% Last modified Feb. 2021

%% Setup environment

% Input checks
if ~exist('task', 'var')
	error('Specify task to be either "yesno" or "yn_threshold"')
end

if ~exist('ft_path', 'var')
	ft_path =  fullfile('fieldtrip20200906', 'fieldtrip-master');
end

if ~exist('data_dir', 'var')
	data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
end

if ~exist('roi', 'var')
	warning('No ROI specified, using default')
    src_dir = fullfile(data_dir, 'source_analysis');
    load(fullfile(src_dir, 'ROI_1f2f_vs_all_noise_in_coi.mat'));
    roi = lf_vox_for_centroid{1};
end

if isempty(roi)
	error('Specify "roi" as valid N x 3 matrix of leadfield voxel coordinates, N > 0')
end


% Environment setup
try
    addpath(genpath('dfi'))
    dfi_startup
catch
    warning('Cannot find dfi folder')
end
addpath(ft_path);
ft_defaults

src_dir = fullfile(data_dir, 'source_analysis');

if strcmp(task, 'yesno')
	an_fold = 'src_fslide_yesno';
elseif strcmp(task, 'yn_threshold')
	an_fold = 'src_fslide_ynt';
end

subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
			'715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N        = length(subjvect);
eegfile  = 'data_preproc2500_to4500.mat';
Fs       = 256;  % Downsample to this frequency 
Nvchan   = size(roi, 1);


%% Load data

% Get leadfield and lcmv_pre (dudd for plotting later)
load(fullfile(src_dir, '701', 'sess2', ...
	'inv_sols_lcmv_c27_100to300ms_600to100msBl_1F2F.mat'), 'leadfield');


%% Project pre-stim data into source space

for isubj = 1:N
    
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
                cfg.toilim   = [-1.2, 0.7];  % matches window in channel fslide analysis
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
    %% 2.) --- SESSION SELECTION ---
    %
    ncomb = length(unique(sessions));
    sessionIDs = unique(sessions);
    nchan = nan(1, length(unique(sessions)));
    
    for isess = 1:ncomb
        
        sess = ['sess', num2str(sessionIDs(isess))];
        
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
        
        % get filter weights (for source space projection)
        load(fullfile(src_dir, partid, sess, ...
			'inv_sols_lcmv_c27_100to300ms_600to100msBl_1F2F.mat'), ...
            'W');
        
        % Get trial IDs
        if length(unique(beh_stim.resp)) == 2
            [id_stim, idm_stim] = dfi_get_trial_indices( beh_stim, 'no' );
        else
            [id_stim, idm_stim] = dfi_get_trial_indices( beh_stim, 'yes' );
        end
        
        % select all trials for pre-stim
        trlsel = [id_stim.v1; id_stim.v2; ...
                  id_stim.av1; id_stim.fus; ...
                  id_stim.fis; id_stim.av2];
        cfg = [];
        cfg.trials = trlsel;
        eeg_stim = ft_selectdata(cfg, eeg_stim);
        beh_stim = beh_stim(trlsel, :);
        
        % Project into source space
        Ntpts = size(eeg_stim.time{1}, 2);
        Ntrls = length(eeg_stim.trialinfo);
        
        src_oi = nan(Nvchan, Ntpts, Ntrls);
        pos_inside = leadfield.pos(leadfield.inside,:);
        for itrl=1:Ntrls
            src = W * eeg_stim.trial{itrl};
            src_oi(:, :, itrl) = src(ismember(pos_inside, roi, 'rows'),:,:);
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
        
        
        %% Compute fslide
        
        s = dfi_eeg_defaults;
        s.chan      = vchan_labels;
        s.psd       = 0;
        s.tfr       = 0;
        s.plv       = 0;
        s.phase     = 0;
        s.fslide    = 1;
        s.filt      = 'no';    % no initial bp filter, we use a tukeyfilter anyway
        s.alpha     = [6, 14];
        s.foi       = s.alpha;
        s.tukeyfilt = 'firls'; % Cohen uses 'firls' with 15% transition,
        % fieldtrip uses 'firws' with half-attenuation at
        % foi(1) and foi(2)
        s.fsplot    = 0;       % Plot processing steps?
        
        cfg         = [];
        cfg.toilim  = [-1.2 0];
        eeg_pre  = ft_redefinetrial(cfg, eeg_src);
        
        % Zero pad post stim period
        s.padwin  = [-1.2 0.7];
        eeg_zeropad = eeg_pre;
        for itrl = 1:numel(eeg_pre.trial)
            eeg_zeropad.time{itrl}  = [eeg_pre.time{itrl}, ...
                -eeg_pre.time{itrl}(end-1:-1:ceil(length(eeg_pre.time{itrl})/2))];
            eeg_zeropad.trial{itrl} = [eeg_pre.trial{itrl}(:,:), ...
                zeros(Nvchan,ceil(length(eeg_pre.time{itrl})/2))];
        end
        
        % fslide
        s.trls = 'all';  
        [~, eeg_fslide]  = dfi_eeg_analysis(s, eeg_zeropad); 
        
        % Save data
        sess_fold = sprintf('session_%i', sessionIDs(isess));
        mkdir(fullfile(subj_dir, an_fold, sess_fold));
        save(fullfile(subj_dir, an_fold, sess_fold, 'fslide_data_trls_all'), ...
			'eeg_fslide', 'beh_stim', '-v7.3');
        
        % Clean up a bit
        close all
        clear eeg_fslide eeg_zeropad eeg_pre
        
    end % sess loop
    
end % subj loop

return

% eof

