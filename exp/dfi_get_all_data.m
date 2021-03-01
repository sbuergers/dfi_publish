%% 0.) --- SETUP ---
clc
close all
clear all

% Name of directory to save things to
an_fold = 'Cecere_et_al_2015';

% save figures?
save_figures = true;

% experiment script folder
addpath(genpath(uigetdir('Select experiment main folder')))

% experiment data folder
data_dir = 'E:\dfi_experiment_data\eeg_data\experiment';
fig_dir  = 'E:\dfi_experiment_figures';

subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N = 20;

task    = 'yesno';
eegfile = 'data_preproc2500_to4500.mat';

for isubj = 1:N

subj_dir = fullfile(data_dir, subjvect{isubj}, task);
    
% get participant, session and session/run information from user
partid    = subjvect{isubj};
rec_type  = 'eyes_open';
resp_lock = 'no';
stim_lock = 'yes';

% figure directories
if save_figures
    if ~exist(fullfile(fig_dir, an_fold, 'PSD'), 'dir'), mkdir(fullfile(fig_dir, an_fold, 'PSD')); end
    if ~exist(fullfile(fig_dir, an_fold, 'PK'),  'dir'), mkdir(fullfile(fig_dir, an_fold, 'PK'));  end
end

% get all data files and append
cnt = 1;
for isess = 1:10
    for irun = 1:10
        temp_dir = fullfile(subj_dir, sprintf('sess%d', isess),  ...
                            sprintf('run%d', irun), rec_type);
        if exist(fullfile(temp_dir, eegfile), 'file')
            fprintf('Loading data from \n%s\n\n', temp_dir)
            load(fullfile(temp_dir, eegfile));
            eeg_cell{cnt} = data_preproc;  clear data_preproc
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

% do we have eye data for all recordings?
eye_data_present = 0;
if exist('eye_behmatch', 'var')
    eye_data_present = 1;
end
if eye_data_present
    for i = 1:numel(beh_cell)
        if isempty(eye_cell{i})
            eye_data_present = 0;
        end
    end
end

% session name
session = 'all_sessions';
supertitle = ['Participant ', partid, ', all sessions', ', all runs'];


























