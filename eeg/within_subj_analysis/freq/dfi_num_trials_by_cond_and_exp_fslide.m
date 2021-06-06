% Show average and std of number of trials for each condition and 
% experiment used in the fslide analyses.
%
% Parent script(s): 
%   dfi_fslide_taking_sessions_into_account_prep.m
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
% Last modified June 2021


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
save_dir = fullfile(data_dir, 'sdt', 'freq_slide');
supp_save_dir = fullfile(data_dir, 'freq_slide');

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

% Which channels do we want to investigate?
channels_of_interest = {'O1', 'O2', 'Oz', 'POz', 'PO4', 'PO8'};


%% Get data

% Yes-no
load(fullfile(data_dir, 'fslide_table_yesno.mat'));
btable_yesno = btable; clear btable

% Delete too fast trials
trls_to_delete = btable_yesno.RT < 0.1;
btable_yesno(trls_to_delete,:) = [];

% Yes-no threshold
load(fullfile(data_dir, 'fslide_table.mat'));
btable_ynt = btable; clear btable

% Delete too fast trials
trls_to_delete = btable_ynt.RT < 0.1;
btable_ynt(trls_to_delete,:) = [];


%% Trial distribution

partvect = unique(btable_ynt.partid);
soavect = unique(btable_yesno.soa);
condvect_1f = [2, 5, 8];
condvect_2f = [3, 6, 9];
[n_1f_yesno, n_2f_yesno, n_1f_ynt, n_2f_ynt] = deal(nan(20, 3));
for isubj = 1:20
    for icond = 1:3
        if icond == 3
            n_1f_yesno(isubj, icond) = sum( ...
                btable_yesno.partid == partvect(isubj) & ...
                btable_yesno.trlid == condvect_1f(icond) & ...
                btable_yesno.soa > 0.042 & btable_yesno.soa < 0.15);
        else
            n_1f_yesno(isubj, icond) = sum( ...
                btable_yesno.partid == partvect(isubj) & ...
                btable_yesno.trlid == condvect_1f(icond));
        end
        n_2f_yesno(isubj, icond) = sum( ...
            btable_yesno.partid == partvect(isubj) & ...
            btable_yesno.trlid == condvect_2f(icond) & ...
            btable_yesno.soa > 0.042 & btable_yesno.soa < 0.15);
        n_1f_ynt(isubj, icond) = sum( ...
            btable_ynt.partid == partvect(isubj) & ...
            btable_ynt.trlid == condvect_1f(icond));
        n_2f_ynt(isubj, icond) = sum( ...
            btable_ynt.partid == partvect(isubj) & ...
            btable_ynt.trlid == condvect_2f(icond));
    end
end

num_sounds = [0; 1; 2];
fprintf('\n\n-----------------------------------------------------------\n\n')
fprintf('       YES-NO  // 1 FLASH\n\n')
disp(table(num_sounds, mean(n_1f_yesno)', std(n_1f_yesno)', std(n_1f_yesno)' ./ sqrt(20),  ...
    'VariableNames', {'num_sounds', 'mean', 'std', 'stderr'}))

fprintf('\n\n-----------------------------------------------------------\n\n')
fprintf('       YES-NO  // 2 FLASHES\n\n')
disp(table(num_sounds, mean(n_2f_yesno)', std(n_2f_yesno)', std(n_2f_yesno)' ./ sqrt(20),  ...
    'VariableNames', {'num_sounds', 'mean', 'std', 'stderr'}))

fprintf('\n\n-----------------------------------------------------------\n\n')
fprintf('       YES-NO THRESHOLD  // 1 FLASH\n\n')
disp(table(num_sounds, mean(n_1f_ynt)', std(n_1f_ynt)', std(n_1f_ynt)' ./ sqrt(20),  ...
    'VariableNames', {'num_sounds', 'mean', 'std', 'stderr'}))

fprintf('\n\n-----------------------------------------------------------\n\n')
fprintf('       YES-NO THRESHOLD  // 2 FLASHES\n\n')
disp(table(num_sounds, mean(n_2f_ynt)', std(n_2f_ynt)', std(n_2f_ynt)' ./ sqrt(20), ...
    'VariableNames', {'num_sounds', 'mean', 'std', 'stderr'}))


% eof 

