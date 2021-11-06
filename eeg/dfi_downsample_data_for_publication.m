% Unfortunately, we were not able to find a solution where we could upload
% the data at the original sampling frequency (1000Hz), and as it was used
% throughout all scripts in this repository. Here we downsample the data
% to 76 Hz, low-pass filter at 28 Hz and cut-out the window between 1200
% ms prior to and 700 ms post-first stimulus.
%
% This step somewhat undermines the replicability of the findings, but it
% should really change them enough to be of huge concern. Note that you
% will need to adapt all scripts in this repository to accept this new data
% format (if this has not been done already).
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
% Last modified Oct. 2021


%% 0.) --- SETUP ---
clc
close all
clear all

% experiment script folder
try
    addpath(genpath('dfi'))
    dfi_startup
catch
    warning('Cannot find dfi folder')
end


% experiment data folder
data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
root_dir_small_data  = 'dfi_small_data';

% add fieldtrip folder to search path
try
    addpath(fullfile('fieldtrip-20160816'))
catch
    warning('Cannot find fieldtrip folder')
end

subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N        = 20;
task     = 'yn_threshold';  % 'yesno', 'yn_threshold'
eegfile  = 'data_preproc2500_to4500.mat';
eegfile_out = 'data_preproc1200to700.mat';
Fs       = 76;  % Downsample to this frequency
lp_fs    = 28;  % Low pass filter at this frequency


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
                
                
                %% 1.) --- DOWNSAMPLE ---
                
                fprintf('Loading data from \n%s\n\n', temp_dir)

                load(fullfile(temp_dir, eegfile));
                load(fullfile(temp_dir, 'artifact_info.mat'));
                                
                % Filter
                cfg              = [];
                cfg.lpfilter     = 'yes';
                cfg.lpfreq       = lp_fs;
                cfg.lpfiltord    = 10;
                cfg.lpfiltdir    = 'twopass';
                data_lowpass     = ft_preprocessing(cfg, data_preproc);

                % Downsample and demean
                cfg            = [];
                cfg.resamplefs = Fs;
                cfg.detrend    = 'no';
                fprintf('\n\nDownsampling to %i Hz...\n\n', cfg.resamplefs);
                stim_downsamp  = ft_resampledata(cfg, data_lowpass);

                % detrend
                cfg          = [];
                cfg.detrend  = 'yes';
                cfg.overlap  = 0;
                data_detrend = ft_preprocessing(cfg, stim_downsamp);
                                
                % select trial window
                cfg          = [];
                cfg.toilim   = [-1.2, 0.7];
                data_cut = ft_redefinetrial(cfg, data_detrend);
                
                % Remove artifacts
                % It seems that simply adjusting the artifact sampling rate
                % does not work for fieldtrip, and generating a sampleinfo
                % after downsampling also does not work very easily. So
                % we will remove artifacts at this stage already. This will
                % potentially change the outcomes of some analyses
                % slightly. In addition, it will break scripts that perform
                % artifact rejection (really all scripts performing
                % preprocessing steps have to be rewritten, sorry about
                % that!).
                
                % delete trials that contain artifacts
                cfg                           = [];
                cfg.artfctdef.muscle.artifact = artifact_samples;
                cfg.artfctdef.reject          = 'complete';
                
                stim_clean = ft_rejectartifact(cfg, data_cut);
        
                % adjust behavioral data accordingly
                oldsampinfo = data_cut.sampleinfo(:,1);
                newsampinfo = stim_clean.sampleinfo(:,1);
                keepTrials  = find( ismember(oldsampinfo, newsampinfo) );

                bdata_stim = bdata(keepTrials, :);
                bdata = bdata_stim;

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
                
                data_preproc = stim_clean;

                % save
                mkdir(fullfile(root_dir_small_data, temp_dir));
                save(fullfile(root_dir_small_data, temp_dir, eegfile_out), ...
                    'data_preproc', 'bdata', '-v7.3');
            end
        end
    end
    
end % loop over participants

% eof

