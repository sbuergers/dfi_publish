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


%% Main loop

% In order to obtain components over which it makes sense to average we
% need to deal with the fact that the orientation in source space (i.e. the
% sign of the ERP) is arbitrary.
% To do so, 
% 1. we will perform a singular value decomposition over voxels and
% consider the first eigenvariate for further analysis.
% 2. next, we will perform another svd over sessions,
% 3. then, we will perform another svd over participants, following the
% same principle. 
% 4. Finally, we obtain SVD eigenvariates over participants
%
% Note that this method, due to sign-flipping, will amplify noise, but
% there is no reason so assume this noise would be ERP shaped. And we are
% interested in visually comparing the source and scalp ERPs.
% TODO ...
% t x nvox
% svd()
% take first eigenvariate
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

eeg_roi = squeeze(nanmean(eeg_1f_src(voxoi, :, :), 3));

img = imread('E:\dfi_experiment_figures\Paper_figures\iAF\src\ortho_slices_at_centroid_of_roi_low_qual.png');
gray = rgb2gray(img);
[U,S,V] = svd(im2double(gray));
r=10;
approx_img = U(:,1:r)*S(1:r,1:r)*V(:,1:r).';
figure; 
subplot(221); image(approx_img)
subplot(222); image(gray)
subplot(223); image(img)

% perform svd over voxels
[U,S,V] = svd(eeg_roi)

X = eeg_roi;
[U,S,V] = svd(X)

figure
plot(S(1,:)); hold on
plot(X(1,:), '-r');
US = U*S;
plot(US(:,1));


% plot
figure
SV = S*V;
plot(SV(1, :)); hold on
plot(X(1, :), '-k', 'linewidth', 2);

% compare to PCA using eigenvalue decomposition
X_norm = X';
X_norm = (X_norm - mean(X_norm))./std(X_norm);
correlMat = corr(X_norm);
[Ve, E] = eig(correlMat);
[coeff, ~, latent] = pca(X_norm);
[sort(diag(E)), sort(testlatent)]

VeX = X'*Ve;
SV = S*V;
USV = U*S*V;

figure;
subplot(221)
title('PC1 and PC2 using eig');
plot(VeX(:,1), 'b'); hold on
plot(VeX(:,2), 'r');
subplot(222)
title('PC2 and PC2 using svd');
plot(SV(1,:), 'b'); hold on
plot(SV(2,:), 'r');
subplot(223)
plot(USV(1,:), 'b'); hold on
plot(USV(2,:), 'r');
subplot(224)
plot(X(1,:), 'b'); hold on
plot(X(2,:), 'r');










load(fullfile(src_dir, '701', 'sess3', 'erps_min600to300ms_1F2F.mat'), ...
    'eeg_1f_trls_avg');
[erps_src, erps_scalp] = deal( nan(N, 10, size(eeg_1f_trls_avg.avg,2)) ); 
choi = ismember(eeg_1f_trls_avg.label, 'O2') | ...
       ismember(eeg_1f_trls_avg.label, 'PO4') | ...
       ismember(eeg_1f_trls_avg.label, 'PO8');

for isubj = 1:N
    
    partid = subjvect{isubj};
    
    fprintf('\n-----------------------------------')
    fprintf('\nParticipant %s\n', partid)
    fprintf('\n         Available sessions \n')
    ls(fullfile(src_dir, sprintf('%s\\sess*',partid)))
    fprintf('-----------------------------------\n')
    
    sessions = ls(fullfile(src_dir, sprintf('%s\\sess*',partid)));
        
    for isess = 1:size(sessions,1)
        
        sess = sessions(isess,1:5);
        
        load(fullfile(src_dir, partid, sess, ...
            'erps_min600to300ms_1F2F.mat'), ...
            'eeg_1f_src_avg', 'eeg_1f_trls_avg');
        
        erps_scalp(isubj, isess, :) = squeeze(nanmean(eeg_1f_trls_avg.avg(choi,:)));
        erps_src(isubj, isess, :) = eeg_1f_src_avg;
        
    end % sess loop

end % subj loop


% Figure of scalp and source evoked response
fh = figure;
time = eeg_1f_trls_avg.time;
subplot(211)
plot(time, squeeze(nanmean(nanmean(erps_scalp,2),1)));
ylabel('Amplitude (uV)')
xlabel('Time (s)')
title('Sensor level')
subplot(212)
plot(time, squeeze(nanmean(nanmean(erps_src,2),1)));
ylabel('Amplitude (a.u.)')
title('Source level')
xlabel('Time (s)')
saveas(fh, fullfile(fig_save_dir, 'ERP_sensor_and_src.svg'))

close all

% enable warnings again
warning('on','all')


% eof


