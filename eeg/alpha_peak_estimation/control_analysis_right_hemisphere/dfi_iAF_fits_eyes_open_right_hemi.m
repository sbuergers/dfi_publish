% Alpha peak frequency estimation for pre-stimulus, task related data
% Yes-no Task
% Use alpha peak frequency fitting method by Corcoran et al. (2018)
%
% Dependencies:
% 	restingIAF (Corcoran et al. toolbox), with modified scripts
% 	meanIAF_sb.m
% 	restingIAF_sb.m
%
% Parent script(s): 
%
% Children script(s): 
%   dfi_bb_analysis_iAF_vs_threshold_allsess_allPFfits_eyes_open.m
%
% Sibling script(s):
%   dfi_iAF_fits_eyes_open_ynt.m
%	dfi_iAF_fits_eyes_closed_w_0padding.m
%
%
% DETAILS
%
% Leverages all posterior channels to obtain a good alpha peak
% frequency estimate in an objective manner.
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

% Name of directory to save things to
an_fold = 'iAF_fits_corcoran_zeropadded';

% save figures?
save_figures = false;

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
task     = 'yesno';
eegfile  = 'data_preproc2500_to4500.mat';

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
tlen   = ceil(0.7*Fs); % window length (corresponds to 0.7 s)
tover  = 0;        % overlap between windows for welch
nfft   = 8*Fs;     % length to zero pad the data to
norm   = true;     % normalize spectra by average power


for isubj = 1:N %[5, 8, 15, 19]
    
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
                
                
                %% 1.) --- PREPROCESSING ---
                
                fprintf('Loading data from \n%s\n\n', temp_dir)

                load(fullfile(temp_dir, eegfile));
                load(fullfile(temp_dir, 'artifact_info.mat'));
                
                % select trial window
                cfg          = [];
                cfg.toilim   = [-1, 0.5];
                data_preproc = ft_redefinetrial(cfg, data_preproc);
                
                fprintf('\n\nPreprocessing stimulus locked data...\n\n');

                % delete trials that contain artifacts
                cfg                           = [];
                cfg.artfctdef.muscle.artifact = artifact_samples;
                cfg.artfctdef.reject          = 'complete';

                stim_clean = ft_rejectartifact(cfg, data_preproc);

                % adjust behavioral data accordingly
                oldsampinfo = data_preproc.sampleinfo(:,1);
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

                eeg_cell_stim{cnt} = stim_reref;
                beh_cell_stim{cnt} = bdata_stim;

                clear stim_downsamp stim_reref stim_preproc stim_seg

                cnt = cnt + 1;
                sessions = [sessions, isess];
                runs     = [runs, irun];
            end
        end
    end
    
    if ~exist('eeg_cell_stim', 'var')
        error('Ups. Could not find EEG data ... Change path so that you can see dfi_experiment_data.')
    end
    
    % session name
    session = 'all_sessions';
    supertitle = ['Participant ', partid, ', all sessions', ', all runs'];
    
    
    %  
    %% 2.) --- SESSION SELECTION (yes-no) ---
    %
    ncomb = length(unique(sessions))+1;
    sessionIDs = unique(sessions);
    nchan = nan(1, length(unique(sessions)));
    
    for isess_combi = 1:ncomb
    
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
        [id_stim, idm_stim] = dfi_get_trial_indices( beh_stim, 'no' );
        
        
        %
        %% 3.) --- ALPHA PEAK FREQUENCY ---
        %
        jx = isess_combi;

        fprintf('\n\nCorcoran iAF estimation:')
        fprintf('Subject: %i, session %i, run %i\n\n', isubj, sessions(jx), runs(jx));
    
        % Look only at pre-stimulus alpha between -0.7 and 0
        cfg          = [];
        cfg.toilim   = [-0.7, 0];
        ps_trls      = ft_redefinetrial(cfg, eeg_stim);

        % select all occipito-parietal sensors
        cfg = [];
        cfg.channel = {'O1', 'O2', 'Oz', 'PO9', 'PO7', 'PO3', 'POz', 'PO4', 'PO8', ...
                       'PO10', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8'}; 
        ps_trls = ft_selectdata(cfg, ps_trls);

        channel_vect = ps_trls.label;

        % treat trials as artifically continous segment (pwelch is going to cut it
        % up again exactly in the same way, and here we use no overlap!)
        lentrl = size(ps_trls.trial{1},2);
        fk_cont_data = nan(numel(ps_trls.label), numel(ps_trls.trial)*lentrl);
        for i = 1:numel(ps_trls.trial)
            fk_cont_data(:, (1+(i-1)*lentrl):(lentrl*i) ) = single(ps_trls.trial{i});
        end

        nchan(jx) = length(ps_trls.label);

        [pSpec(isubj, jx).sums, pSpec(isubj, jx).chans, f] = restingIAF_sb(fk_cont_data, nchan(jx), cmin, fRange, Fs, w, Fw, k, mpow, mdiff, tlen, tover, nfft, norm);

    end % session combinations
    
    % weighted average of mean IAF estimates across (j-th) recordings
    if length(unique(sessions)) > 1
        % i.) all runs
        [muPaf(isubj,:), muCog(isubj,:), muPow(isubj,:), ...
         muPow_osc(isubj,:), muPow_pk(isubj,:), muPow_pk_osc(isubj,:)] = meanIAF_sb([pSpec(isubj, 1:length(unique(sessions))).sums], nchan(1:length(unique(sessions))));
    end
    
    clear eeg_cell_stim beh_cell_stim
    
end % loop over participants


figure
plot(muPaf, 'b'); hold on
title('Indivisual alpha peak frequencies (Corcoran method)')
xlabel('Subject ID'); set(gca, 'xtick', 1:20);
ylabel('iAF (Hz)'); grid on; xlim([1 20])


% Single run iAF and periodograms for all subjects (log scale)
selected_channel = 'POz'; %'O2'
kx = strcmp(channel_vect, selected_channel);
figure('color', 'w', 'position', [50 50 1300 900]);
for isubj = 1:N
    subplot(4,5,isubj)
    pChans = pSpec(isubj, 1).chans(kx);
    plot(f, log10(pChans.pxx)); hold on
    plot(f, log10(pChans.d0), 'm'); yl = ylim;
    plot(f, pChans.minPow, 'r'); xlim([3 18]);%xlim(fRange)
    plot(f(f == pChans.peaks), log10(pChans.pxx(f == pChans.peaks)), '.g', 'markersize', 15); hold on
    line([pChans.inf1 pChans.inf1], [yl(1) yl(2)], 'color', 'k');
    line([pChans.inf2 pChans.inf2], [yl(1) yl(2)], 'color', 'k');
    if isubj > 16, xlabel('Frequency'); end
    ylabel('Normalized power'); title(sprintf('subj %i', isubj));
    plotspecs
end
suptitle(sprintf('Channel %s', selected_channel))



% Single run iAF and periodograms for all subjects
selected_channel = 'POz'; %'O2'
kx = strcmp(channel_vect, selected_channel);
figure('color', 'w', 'position', [50 50 1300 900]);
for isubj = 1:N
    subplot(4,5,isubj)
    pChans = pSpec(isubj, 1).chans(kx);
    plot(f, pChans.pxx); hold on
    plot(f, pChans.d0, 'm');
    plot(f, 10.^(pChans.minPow), 'r');
    plot(f(f == pChans.peaks), pChans.d0(f == pChans.peaks), '.g', 'markersize', 15); hold on
    yl = ylim; xlim([3 18]);%xlim(fRange)
    line([pChans.inf1 pChans.inf1], [yl(1) yl(2)], 'color', 'k');
    line([pChans.inf2 pChans.inf2], [yl(1) yl(2)], 'color', 'k');
    if isubj > 16, xlabel('Frequency'); end
    ylabel('Normalized power'); title(sprintf('subj %i', isubj));
    plotspecs
end
suptitle(sprintf('Channel %s', selected_channel))
    


%% 4.) --- SAVE DATA ---

if ~exist(fullfile(fig_dir, an_fold), 'dir')
    mkdir(fullfile(fig_dir, an_fold));
end

fname = 'eyes_open_pkinfo';
save(fullfile(fig_dir, an_fold, fname), ...
    'muPaf*', 'muPow*', 'pSpec', 'cmin', 'channel_vect', ...
    'fRange', 'w', 'Fs', 'Fw', 'k');


% turn matlab warnings back on
warning('on','all')

close all


% eof

