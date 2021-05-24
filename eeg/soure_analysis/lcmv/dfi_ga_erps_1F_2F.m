% Plot grand average scalp and source ERPs for 1F and 2F trials.
%
% Parent script(s): 
%   dfi_erps_1F_2F.m
%
% Children script(s): 
%   None
%
% Sibling script(s):
%   None
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
fig_dir = fullfile('dfi_experiment_figures', 'Paper_figures');
fig_save_dir = fullfile(fig_dir, 'iAF', 'src'); 

subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N        = length(subjvect);
Fs       = 64;  % Downsample to this frequency
tasks    = {'yesno', 'yn_threshold'};


% In order to obtain components over which it makes sense to average we
% need to deal with the fact that the orientation in source space (i.e. the
% sign of the ERP) is arbitrary.
% To do so, 
% 1. we will perform a singular value decomposition over voxels and
% consider the first eigenvariate for further analysis.
% 2. next, we will perform another svd over sessions, following the
% same principle. 
% 3. Finally, we will flip the orientation of the first eigenvariate of
% each participant if the MSE of it vs the scalp ERP is smaler than the MSE
% of scalp ERP and unflipped eigenvariate (this biases the results somewhat
% toward `fitting` the ERP's shape).

% Get voxels of interest
load(fullfile(src_dir,  '701', 'sess3', ...
    'inv_sols_lcmv_c27_100to300ms_600to100msBl_1F2F.mat'), ...
    'W', 'leadfield');
load(fullfile(src_dir, 'ROI_1f2f_vs_all_noise_in_coi.mat'), ...
    'lf_vox_for_centroid', 'centroids');
voxoi = nan(size(lf_vox_for_centroid{1}, 1), 1);
for ivox = 1:size(lf_vox_for_centroid{1}, 1)
    voxoi(ivox) = find(ismember(leadfield.pos, ...
        lf_vox_for_centroid{1}(ivox, :),'rows'));
end
Nvox = length(voxoi);

% Initialize variables
load(fullfile(src_dir, '701', 'sess3', 'erps_min600to300ms_1F2F.mat'), ...
    'eeg_1f_trls_avg');
Nt = size(eeg_1f_trls_avg.avg,2);
[erps_src, erps_scalp] = deal( nan(N, 10, Nt) ); 
choi = ismember(eeg_1f_trls_avg.label, 'O2') | ...
       ismember(eeg_1f_trls_avg.label, 'PO4') | ...
       ismember(eeg_1f_trls_avg.label, 'PO8');
eigenvar_subj = nan(N, Nt);
singval_subj = nan(N, 1);

% Main loop
for isubj = 1:N
    
    partid = subjvect{isubj};
    
    fprintf('\n-----------------------------------')
    fprintf('\nParticipant %s\n', partid)
    fprintf('\n         Available sessions \n')
    ls(fullfile(src_dir, sprintf('%s\\sess*',partid)))
    fprintf('-----------------------------------\n')
    
    sessions = ls(fullfile(src_dir, sprintf('%s\\sess*',partid)));
    Nsess = size(sessions,1);
        
    eigenvar_sess = nan(Nsess, Nt); 
    singval_sess = nan(Nsess, 1);
    for isess = 1:Nsess
        
        sess = sessions(isess,1:5);
        
        load(fullfile(src_dir, partid, sess, ...
            'erps_min600to300ms_1F2F.mat'), ...
            'eeg_1f_src', 'eeg_1f_trls_avg');
        
        erps_scalp(isubj, isess, :) = mean(eeg_1f_trls_avg.avg(choi,:));
        
        eeg_roi = squeeze(nanmean(eeg_1f_src(voxoi, :, :), 3));

        % SVD 1: Over voxels
        X = eeg_roi;
        [U,S,V] = svd(X);
        eigenvar_sess(isess,:) = V(:,1);
        singval_sess(isess) = S(1);
        
%         % Visually compare first eigenvariate and first voxel's signal
%         subplot(3,3,isess)
%         plot(V(:,isess)./max(V(:,isess))); hold on
%         plot(X(isess,:)./max(X(isess,:)));
        
    end % sess loop
    
    % SVD 2: Over sessions
    X = eigenvar_sess;
    [U,S,V] = svd(X);
    eigenvar_subj(isubj,:) = V(:,1);
    singval_subj(isubj) = S(1);
    
%     % Plot
%     subplot(3,3,6)
%     plot(V(:,1)./max(V(:,1))); hold on
%     plot(X(1,:)./max(X(1,:)));
    
end % subj loop

% scalp ERP for each subj
erps_subj = squeeze(nanmean(erps_scalp, 2));

% SVD 3: Over participants
X = eigenvar_subj;
[U,S,V] = svd(X);

% Compute mean squared error between scalp ERPs and
% each participant's first eigenvariate over sessions either flipped
% or not. Choose the one that minimizes the MSE. This way we should get
% consistent orientations over participants and should be able to average.
% Of course, this somewhat biases the results toward fitting with the ERP.
eigenvar_subj_adj = eigenvar_subj;
erp_twin = time >= -0.7;
for isubj = 1:N
    mse_pos = mean((erps_subj(isubj,erp_twin) - eigenvar_subj(isubj,erp_twin)).^2);
    mse_neg = mean((erps_subj(isubj,erp_twin) + eigenvar_subj(isubj,erp_twin)).^2);
    if mse_neg < mse_pos
        eigenvar_subj_adj(isubj,:) = -eigenvar_subj(isubj,:);
    end
end

% Plot individual eigenvariates for participants vs their scalp ERPs
figure
for isubj = 1:20
    subplot(4,5,isubj)
    plot(erps_subj(isubj,:) ./ max(abs(erps_subj(isubj,:)))); hold on
    plot(eigenvar_subj(isubj,:) ./ max(abs(eigenvar_subj(isubj,:))));
    plot(eigenvar_subj_adj(isubj,:) ./ max(abs(eigenvar_subj_adj(isubj,:))));
end

% Compute grand average over eigenvariates
eigenvar_ga = nanmean(eigenvar_subj_adj);

% Figure of scalp and source evoked response
fh = figure;
time = eeg_1f_trls_avg.time;
subplot(211)
plot(time, squeeze(nanmean(nanmean(erps_scalp,2),1)));
ylabel('Amplitude (uV)')
xlabel('Time (s)')
title('Sensor level')
subplot(212)
plot(time, eigenvar_ga);
ylabel('Amplitude (a.u.)')
title('Source level')
xlabel('Time (s)')
fh.Renderer = 'painters'; 
mkdir(fig_save_dir)
saveas(fh, fullfile(fig_save_dir, 'ERP_sensor_and_src.svg'))

close all

% enable warnings again
warning('on','all')


% eof


