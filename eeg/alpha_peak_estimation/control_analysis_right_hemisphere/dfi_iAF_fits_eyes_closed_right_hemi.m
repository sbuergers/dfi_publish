% Alpha peak frequency estimation for eyes-closed resting data
% Yes-no threshold task
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
%   dfi_bb_analysis_iAF_vs_threshold_allsess_eyes_closed.m
%
% Sibling script(s):
%   dfi_iAF_fits_eyes_open.m
%	dfi_iAF_fits_eyes_open_ynt.m
%
%
% DETAILS
%
% Leverages only right hemispheric posterior channels to obtain a good 
% alpha peak frequency estimate in an objective manner.
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

% turn off matlab warnings
warning('off','all')

% Name of directory to save things to
an_fold = 'iAF_fits_corcoran_zeropadded';

% save figures?
save_figures = true;

% experiment data folder
data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
fig_dir  = fullfile('dfi_experiment_figures');

% experiment script folder
try
    addpath(genpath('dfi'))
    dfi_startup
catch
    warning('Cannot find dfi folder')
end

% run startup function
dfi_startup

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
eegfile  = 'eeg_data_eyes_closed.mat';

% Initialize variables for iAF estimation
Fs   = 256; % Downsample to this frequency
[pSpec]    = deal(struct('chans', [], 'sums', []));
muPaf      = nan(N, 1);
muPaf_pre  = nan(N, 1);
muPaf_post = nan(N, 1);
muCog      = nan(N, 1);
muCog_pre  = nan(N, 1);
muCog_post = nan(N, 1);
muPow_abs      = nan(N, 1);
muPow_abs_pre  = nan(N, 1);
muPow_abs_post = nan(N, 1);
muPow_osc      = nan(N, 1);
muPow_osc_pre  = nan(N, 1);
muPow_osc_post = nan(N, 1);
muPkPow_abs      = nan(N, 1);
muPkPow_abs_pre  = nan(N, 1);
muPkPow_abs_post = nan(N, 1);
muPkPow_osc      = nan(N, 1);
muPkPow_osc_pre  = nan(N, 1);
muPkPow_osc_post = nan(N, 1);
sessionmat = nan(N, 11);

% Peak fitting parameter settings
cmin   = 3;        % min # of channels required for cross-channel average
fRange = [3 30];   % frequency range going into fitting process
w      = [6 14];   % alpha peak search window
Fw     = 11;       % Window of polinomial filter 
k      = 5;        % Order of polinomial filter (SGT)
mpow   = 1;        % error bound used to determine threshold differentiating 
                   % substantive peaks from background spectral noise 
mdiff  = 0.2;      % minimum difference between peaks
tlen   = 4*Fs;     % window length (corresponds to 2 seconds)
tover  = 0.5*tlen; % overlap between windows for welch
nfft   = 8*Fs;     % will use next power of 2 by default
norm   = true;     % normalize spectra by average power

for isubj = 1:N
    
    clear eeg_cell eeg_cell_stim

    subj_dir = fullfile(data_dir, subjvect{isubj}, task);

    % get participant, session and session/run information from user
    partid    = subjvect{isubj};
    rec_type  = 'eyes_closed';
    resp_lock = 'no';
    stim_lock = 'yes';

    % figure directories
    if save_figures
        if ~exist(fullfile(fig_dir, an_fold, 'PSD'), 'dir'), mkdir(fullfile(fig_dir, an_fold, 'PSD')); end
        if ~exist(fullfile(fig_dir, an_fold, 'PK'),  'dir'), mkdir(fullfile(fig_dir, an_fold, 'PK'));  end
    end

    % get all data files and append
    % 1.) yes-no sessions
    sessions = [];
    runs     = [];
    for isess = 1:10
        for irun = 1:10
            temp_dir = fullfile(subj_dir, sprintf('sess%d', isess),  ...
                                sprintf('run%d', irun), rec_type);
            if exist(fullfile(temp_dir, eegfile), 'file')
                fprintf('Loading data from \n%s\n\n', temp_dir)
                load(fullfile(temp_dir, eegfile));
                eeg_cell{isess}{irun} = data_clean;  clear data_clean
                sessions = [sessions, isess];
                runs = [runs, irun];
            end
        end
    end
    
    sessionmat(isubj,1:length(sessions)) = sessions;
    
    % 2.) yes-no threshold sessions
    sessions2 = [];
    runs2     = [];
    yn_th_dir = fullfile(data_dir, subjvect{isubj}, 'yn_threshold');
    for isess = 1:10
        for irun = 1:10
            temp_dir = fullfile(yn_th_dir, sprintf('sess%d', isess),  ...
                                sprintf('run%d', irun), rec_type);
            if exist(fullfile(temp_dir, eegfile), 'file')
                fprintf('Loading data from \n%s\n\n', temp_dir)
                load(fullfile(temp_dir, eegfile));
                eeg_cell2{isess}{irun} = data_clean;  clear data_clean
                sessions2 = [sessions2, isess];
                runs2 = [runs2, irun];
            end
        end
    end
    
    sessionmat(isubj,sum((~isnan(sessionmat(isubj,:))),2)+1: ...
                     sum((~isnan(sessionmat(isubj,:))),2)+length(sessions2)) = sessions2;
    
    


    %% 1.) --- PREPROCESSING ---

    for ifile = 1:length(sessions)

        try
            fprintf('\n\nPreprocessing stimulus locked data...\n\n');

            % Downsample and demean
            cfg            = [];
            cfg.resamplefs = Fs;
            cfg.detrend    = 'no';
            fprintf('\n\nDownsampling to %i Hz...\n\n', cfg.resamplefs);
            stim_downsamp  = ft_resampledata(cfg, eeg_cell{sessions(ifile)}{runs(ifile)});

            % Rereference
            fprintf('\n\nRereferencing...\n\n');
            cfg            = [];
            cfg.channel    = {'all'; '-blink'; '-pupil'; '-gazeX'; '-gazeY'};
            cfg.reref      = 'yes';
            cfg.refchannel = 'all';
            stim_reref     = ft_preprocessing(cfg, stim_downsamp);

            eeg_cell_stim{sessions(ifile)}{runs(ifile)} = stim_reref;

            clear stim_blcorr stim_downsamp stim_reref stim_clean ...
                eye_stim stim_preproc stim_filt
        catch
        end
            
    end % file loop

    clear eeg_cell 
    
    
    
    for ifile = 1:length(sessions2)

        try
            fprintf('\n\nPreprocessing stimulus locked data...\n\n');

            % Downsample and demean
            cfg            = [];
            cfg.resamplefs = Fs;
            cfg.detrend    = 'no';
            fprintf('\n\nDownsampling to %i Hz...\n\n', cfg.resamplefs);
            stim_downsamp  = ft_resampledata(cfg, eeg_cell2{sessions2(ifile)}{runs2(ifile)});

            % Rereference
            fprintf('\n\nRereferencing...\n\n');
            cfg            = [];
            cfg.channel    = {'all'; '-blink'; '-pupil'; '-gazeX'; '-gazeY'};
            cfg.reref      = 'yes';
            cfg.refchannel = 'all';
            stim_reref     = ft_preprocessing(cfg, stim_downsamp);

            eeg_cell_stim{sessions2(ifile)}{runs2(ifile)} = stim_reref;

            clear stim_blcorr stim_downsamp stim_reref stim_clean ...
                eye_stim stim_preproc stim_filt
        catch
        end
            
    end % file loop

    clear eeg_cell2 
    
    sessions = [sessions, sessions2];
    runs = [runs, runs2];


    
    %% 3.) --- ALPHA PEAK FREQUENCY (Corcoran, 2017) ---
    
    nchan = nan(1, length(sessions));
    
    for jx = 1:length(sessions)
        
        fprintf('\n\nCorcoran iAF estimation:')
        fprintf('Subject: %i, session %i, run %i\n\n', isubj, sessions(jx), runs(jx));
        
        eeg_stim = eeg_cell_stim{sessions(jx)}{runs(jx)};
        
        % Cut data into 2 s pieces (without nans)
        cfg              = [];
        cfg.length       = 2;
        cfg.minlength    = 2;
        cfg.overlap      = 0;
        ps_trls          = ft_redefinetrial(cfg, eeg_stim);
        
        % select all right-hemispheric occipito-parietal sensors
        cfg = [];
        cfg.channel = {'O2', 'Oz', 'POz', 'PO4', 'PO8', ...
                       'PO10', 'Pz', 'P2', 'P4', 'P6', 'P8'}; 
        ps_trls = ft_selectdata(cfg, ps_trls);
        
        channel_vect = ps_trls.label;
        
        % treat trials as artifically continous segment (pwelch is going to cut it
        % up again exactly in the same way, except that it also uses 0.5s overlap.
        lentrl = size(ps_trls.trial{1},2);
        fk_cont_data = nan(numel(ps_trls.label), numel(ps_trls.trial)*lentrl);
        for i = 1:numel(ps_trls.trial)
            fk_cont_data(:, (1+(i-1)*lentrl):(lentrl*i) ) = single(ps_trls.trial{i});
        end

        nchan(jx) = length(ps_trls.label);
        
        [pSpec(isubj, jx).sums, pSpec(isubj, jx).chans, f] = restingIAF_sb(fk_cont_data, nchan(jx), cmin, fRange, Fs, w, Fw, k, mpow, mdiff, tlen, tover, nfft, norm);
        
    end
    
    % weighted average of mean IAF estimates across (j-th) recordings
    if length(sessions) > 1
        % i.) all runs
        [muPaf(isubj,:), muCog(isubj,:), muPow_abs(isubj,:), ...
         muPow_osc(isubj,:), muPkPow_abs(isubj,:), muPkPow_osc(isubj,:)] = ...
                                                    meanIAF_sb([pSpec(isubj,:).sums], nchan);
    end
    
    if length(sessions) > 1 && length(unique(runs)) > 1
        % ii.) only pre-stimulus runs
        [muPaf_pre(isubj,:), muCog_pre(isubj,:), muPow_abs_pre(isubj,:), ...
         muPow_osc_pre(isubj,:), muPkPow_abs_pre(isubj,:), muPkPow_osc_pre(isubj,:)] = ...
                                                    meanIAF_sb([pSpec(isubj, runs==1).sums], nchan(runs==1));
        % iii.) only post-stimulus runs
        [muPaf_post(isubj,:), muCog_post(isubj,:), muPow_abs_post(isubj,:), ...
         muPow_osc_post(isubj,:), muPkPow_abs_post(isubj,:), muPkPow_osc_post(isubj,:)] = ...
                                                    meanIAF_sb([pSpec(isubj, runs==2).sums], nchan(runs==2));
    end          
    
    clear eeg_cell_stim beh_cell_stim

end % loop over participants

[muPaf, muPaf_pre, muPaf_post]
figure
plot(muPaf_pre, 'b'); hold on
plot(muPaf_post, 'r'); plotspecs
title('Indivisual alpha peak frequencies (Corcoran method)')
xlabel('Subject ID'); set(gca, 'xtick', 1:20);
ylabel('iAF (Hz)'); grid on; xlim([1 20])


% Single run iAF and periodograms for all subjects (log scale)
selected_channel = 'PO4'; %'O2'
kx = strcmp(channel_vect, selected_channel);
figure('color', 'w', 'position', [50 50 1300 900]);
f = fRange(1):0.125:fRange(end);
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
selected_channel = 'PO4'; %'O2'
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




% Single run iAF and periodograms for all subjects (log scale)
% All channels!
channels_to_show = {'O2', 'PO4', 'PO8'};
selected_channel = channels_to_show{2}; 
kx = strcmp(channel_vect, selected_channel);
figure('color', 'w', 'position', [50 50 900 600]);
ha = tight_subplot(10, 20,[0.005 0.005],[0.01],[0.01]);
f = fRange(1):0.125:fRange(end);
subidx = 1:20;
recidx = 1:10;
for isubj = 1:N
    for irec = 1:10
    %for ichan = 1:length(channels_to_show)
        %axes(ha(10*(isubj-1) + irec));
        axes(ha(isubj + 20*(irec-1)));
        try
        pChans = pSpec(isubj, irec).chans(kx);
        plot(f, log10(pChans.pxx)); hold on
        plot(f, log10(pChans.d0), 'm'); yl = ylim;
        plot(f, pChans.minPow, 'r');
        plot(f(f == pChans.peaks), log10(pChans.pxx(f == pChans.peaks)), '.k', 'markersize', 12); hold on
        line([pChans.inf1 pChans.inf1], [yl(1) yl(2)], 'color', 'k');
        line([pChans.inf2 pChans.inf2], [yl(1) yl(2)], 'color', 'k');
        ylim([-1 2])
        xlim([5 15])
        set(gca, 'xticklabel', [])
        set(gca, 'yticklabel', [])
        catch
        end
    end
end
suptitle(sprintf('Channel %s', selected_channel))




%% 4.) --- SAVE DATA ---

if ~exist(fullfile(fig_dir, an_fold), 'dir')
    mkdir(fullfile(fig_dir, an_fold));
end

fname = 'eyes_closed_pkinfo_yn_plus_ynt';
save(fullfile(fig_dir, an_fold, fname), ...
    'f', 'muPaf*', 'muPow*', 'muPkPow*', 'pSpec', ...
    'cmin', 'channel_vect', 'fRange', 'w', 'Fs', 'Fw', 'k');


% turn matlab warnings back on
warning('on','all')

close all


% // eof

















