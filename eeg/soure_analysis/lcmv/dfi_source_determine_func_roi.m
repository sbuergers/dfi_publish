% Determine functional ROI from leadfield voxels
%
% Parent script(s): 
%   dfi_source_proj_erp_group_lvl.m
%
% Children script(s): 
%   dfi_source_fslide_prep.m
%
% Sibling script(s):
%   None
%
%
% DETAILS
%
% The region with the strongest effect tends to also have relatively large
% variability over subjects. Nevertheless, the grand average contrast map 
% seems to correspond better to what we might expect than the tmaps, 
% possibly, because small, but consistent effects over subjects yield
% larger t-values contralateral to stimulus presentation. 
%
% To get a good estimate in the right hemisphere (contralteral to stim.
% presentation), we will therefore pre-select all voxels with significant
% t-values at the 0.05 level (SELECTION 1), and then select the largest grand
% average contrast values within that region (SELECTION 3). In addition,
% we will constrain ourselves to the anatomical ROI specified above (i.e.
% only visual areas in the right hemisphere (SELECTION 2). 
%
% To recap:
% 1.) Constrain ourselves to significant voxels at the 0.05 alpha level
% 2.) Constrain ourselves to anatomical visual ROI in right hemisphere
% 3.) Find largest grand average contrast values within that pre-selection
%
% Use the following contrast (normalized by voxel average)
% source_effect_mean = (mean(var(post_src, 0, 2),3) - ...
%                       mean(var(pre_src, 0, 2),3)) ./ ...
%                       mean(var(pre_src, 0, 2),3);
% source_effect_mean = source_effect_mean - mean(source_effect_mean);
% 
% NOTES
% For this source analysis we use a newer fieldtrip version than for
% the rest of the analyses (quite a bit has changed in fieldtrip 
% regarding sourceanalysis over the years so this simply seems 
% least likely to lead to problems).
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



%% Setup environment

restoredefaultpath; clc; close all; clear all
warning('off','all')

% Add paths
try
    addpath(genpath('dfi'))
    dfi_startup
catch
    warning('Cannot find dfi folder')
end
ft_path = fullfile('fieldtrip20200906', 'fieldtrip-master'); 
addpath(ft_path);
ft_defaults
addpath(fullfile('spm12', 'spm12'))  % for atlas in sourceplots

data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
fig_dir = fullfile('dfi_experiment_figures', 'Paper_figures');
src_dir = fullfile(data_dir, 'source_analysis');
fig_save_dir = fullfile(fig_dir, 'iAF', 'src'); 
mkdir(fig_save_dir);

% Colin27 anatomical
mri = ft_read_mri(fullfile(ft_path, 'template', 'headmodel', 'standard_mri.mat'));


%% Load data
load(fullfile(src_dir, 'fun_localizer_1F2F_vs_all_noise_in_coi.mat'));

% Compute effect with second pipeline: 
% Source space contrast 2:
% pre_src = reshape(W * cell2mat(eeg_pre.trial), [Nvox, Ntpts, Ntrls]);
% post_src = reshape(W * cell2mat(eeg_post.trial), [Nvox, Ntpts, Ntrls]);
%  
% source_effect = mean((var(post_src,0,2)-var(pre_src,0,2))./var(pre_src,0,2),3);
% source_effect_mean = (mean(var(post_src, 0, 2),3) - ...
%                       mean(var(pre_src, 0, 2),3)) ./ ...
%                       mean(var(pre_src, 0, 2),3);

n1 = sum(leadfield.inside);

% Convert contrast cell arrays to matrices
contrast2_mat = reshape(cell2mat(contrast2), [n1, 20]);
contrast2_mat_norm_sub = reshape(cell2mat(contrast2_norm_sub), [n1, 20]);

contrast2_mat_mean = reshape(cell2mat(contrast2_mean), [n1, 20]);
contrast2_mat_mean_norm_sub = reshape(cell2mat(contrast2_mean_norm_sub), [n1, 20]);

% Compute grand average contrasts
GA_contrast2 = mean(reshape(cell2mat(contrast2), [n1, 20]), 2);
GA_contrast2_norm_sub = mean(reshape(cell2mat(contrast2_norm_sub), [n1, 20]), 2);

GA_contrast2_mean = mean(reshape(cell2mat(contrast2_mean), [n1, 20]), 2);
GA_contrast2_mean_norm_sub = mean(reshape(cell2mat(contrast2_mean_norm_sub), [n1, 20]), 2);


% Compute t-values
[tmap_ratio2, tmap_ratio2_norm_sub, ...
    tmap_ratio2_mean, tmap_ratio2_mean_norm_sub] = deal(nan(size(GA_contrast2)));
disp('computing t-values over contrasts..... ')
for ivox = 1:n1
    
    [~, p, ~, stats] = ttest(contrast2_mat(ivox, :));
    tmap_ratio2(ivox) = stats.tstat;
    [~, p, ~, stats] = ttest(contrast2_mat_norm_sub(ivox, :));
    tmap_ratio2_norm_sub(ivox) = stats.tstat;
    
    [~, p, ~, stats] = ttest(contrast2_mat_mean(ivox, :));
    tmap_ratio2_mean(ivox) = stats.tstat;
    [~, p, ~, stats] = ttest(contrast2_mat_mean_norm_sub(ivox, :));
    tmap_ratio2_mean_norm_sub(ivox) = stats.tstat;
end
disp('..... done.')



%% Determine ROI

% Load atlas from Wang et al. (2015)
% Wang, Liang, Ryan E. B. Mruczek, Michael J. Arcaro, and Sabine Kastner. 
% “Probabilistic Maps of Visual Topography in Human Cortex.” Cerebral Cortex 
% (New York, N.Y.: 1991) 25, no. 10 (October 2015): 3911–31. 
% https://doi.org/10.1093/cercor/bhu277.

load(fullfile(ft_path, 'template', 'atlas', 'vtpm', 'vtpm.mat'))

figure;
subplot(2,2,1); imagesc(vtpm.tissue(:,:,70))
subplot(2,2,2); imagesc(vtpm.tissue(:,:,80))
subplot(2,2,3); imagesc(vtpm.tissue(:,:,90))
subplot(2,2,4); imagesc(vtpm.tissue(:,:,100))

disp(vtpm)

% Note that the dimensions of the atlas and the colin27 mri do not quite
% correspond, yet both should in principle be in MNI space already.
% Colin 27: [181 217 181]          vtpm: [182 218 182]

% Interpolate to same space
cfg = [];
cfg.parameter = 'tissue';
mri = ft_sourceinterpolate(cfg, vtpm, mri);
mri.tissuelabel = vtpm.tissuelabel;

% Let's use the whole visual system (without multisensory areas).
% Order of labels should be the same as in vtpm.tissuelabel!
visual_ROI_labels = {'right_V1v', 'right_V1d', 'right_V2v', 'right_V2d', ...
                     'right_V3v', 'right_V3d', 'right_hV4', ...
                     'right_VO1', 'right_VO2', 'right_PHC1', 'right_PHC2', ...
                     'right_MST', 'right_hMT', 'right_LO2', 'right_LO1', ...
                     'right_V3b', 'right_V3a'};
visual_ROI = find(ismember(vtpm.tissuelabel, visual_ROI_labels));


% We will use the following procedure for defining our ROI:
% The region with the strongest effect tends to also have relatively large
% variability over subjects. Nevertheless, the grand average contrast map 
% seems to correspond better to what we might expect than the tmaps, 
% possibly, because small, but consistent effects over subjects yield
% larger t-values contralateral to stimulus presentation. 
%
% To get a good estimate in the right hemisphere (contralteral to stim.
% presentation), we will therefore pre-select all voxels with significant
% t-values at the 0.05 level (SELECTION 1), and then select the largest grand
% average contrast values within that region (SELECTION 3). In addition,
% we will constrain ourselves to the anatomical ROI specified above (i.e.
% only visual areas in the right hemisphere (SELECTION 2). 
%
% To recap:
% 1.) Constrain ourselves to significant voxels at the 0.05 alpha level
% 2.) Constrain ourselves to anatomical visual ROI in right hemisphere
% 3.) Find largest grand average contrast values within that pre-selection
%
% Use the following contrast (normalized by voxel average)
% source_effect_mean = (mean(var(post_src, 0, 2),3) - ...
%                       mean(var(pre_src, 0, 2),3)) ./ ...
%                       mean(var(pre_src, 0, 2),3);
% source_effect_mean = source_effect_mean - mean(source_effect_mean);


% 1.) Pre-selection: t-values > 2.093 [2-tailed 0.05 alpha threshold at N=20]
sig_t = tmap_ratio2_mean_norm_sub > 2.093;

dummy = lcmv_all;
dummy.effect = nan(size(leadfield.inside));
dummy.effect(leadfield.inside) = tmap_ratio2_mean_norm_sub;
dummy.mask = leadfield.inside;
dummy.mask(leadfield.inside) = sig_t;

cfg = [];
cfg.parameter = 'effect';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

cfg = [];
cfg.parameter = 'mask';
[interp2] = ft_sourceinterpolate(cfg, dummy, mri);
interp_mask_sig_t = interp2.mask;  clear interp2

interp.mask = interp_mask_sig_t;
cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'effect';
cfg.maskparameter = 'mask';
cfg.funcolorlim = [0 max(dummy.effect)];
ft_sourceplot(cfg, interp);
colormap('hot')
colormap('viridis')


% Combined pre-selection GA contrast (norm sub)
% (1. t-values > 2.1; 
%  2. in right visual ROI)
dummy = lcmv_all;
dummy.effect = nan(size(leadfield.inside));
dummy.effect(leadfield.inside) = mean(contrast2_mat_mean_norm_sub, 2);

cfg = [];
cfg.parameter = 'effect';
[interp] = ft_sourceinterpolate(cfg, dummy, mri);

interp.mask = reshape(interp_mask_sig_t, [mri.dim]) & ...
              ismember(mri.tissue, visual_ROI);
cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'effect';
cfg.maskparameter = 'mask';
cfg.funcolorlim = [min(dummy.effect) max(dummy.effect)];
ft_sourceplot(cfg, interp);
colormap('hot')
colormap('viridis')
saveas(gcf, fullfile(fig_save_dir, 'right_hemi_contrast_in_atlas_slices.png'))


% Determine the 100*x% top percentile contrast values within the ROI
x = 0.05;
mask = ismember(mri.tissue, visual_ROI);
mask_long = interp_mask_sig_t & ...
            reshape(mask, [size(mask,1)*size(mask,2)*size(mask,3), 1]);

effect_subset = interp.effect(mask_long);
effect_subset = sort(effect_subset(~isnan(effect_subset)));
cutoff = effect_subset(round(length(effect_subset)*(1-x)));

sum(effect_subset > cutoff) / length(effect_subset)  % should be x


% plot grand average effect inside ROI mask with x% largest power 
interp.mask_effect = interp.effect > cutoff & mask_long;
cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'effect';
cfg.maskparameter = 'mask_effect';
cfg.funcolorlim = [min(dummy.effect) max(dummy.effect)];
ft_sourceplot(cfg, interp);
colormap('hot')
colormap('viridis')

% plot 3 orthogonal slices
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'effect';
cfg.maskparameter = 'mask_effect';
cfg.funcolorlim   = [min(dummy.effect) max(dummy.effect)];
cfg.atlas = vtpm;
ft_sourceplot(cfg, interp);
colormap('viridis')
saveas(gcf, fullfile(fig_save_dir, 'roi_ortho_slices.png'))


% Which anatomical regions are included here?
mri_tissue_long = reshape(mri.tissue, [size(mask,1)*size(mask,2)*size(mask,3), 1]);
roi_included = zeros(size(visual_ROI_labels,1),1);
for i = 1:length(visual_ROI_labels)
    roi_included(i) = sum(mri_tissue_long == i & ...  % is in ROI
                          interp.effect > cutoff & mask_long);  % is in selection
end
disp('Which anatomical regions are included in our selection?')
for i = 1:length(visual_ROI_labels)
    fprintf('%s has %i voxels included\n', ...
        visual_ROI_labels{i}, roi_included(i));
end


% It looks like we have three distinct clusters at 10%, but only 1 at 5%:
%
% 10%
% LO1 and LO2
% hV4 and VO1
% V3b
%
% 5%
% hV4 and VO1

% Determine number of clusters by adjacency
CC = bwconncomp(reshape(interp.mask_effect, [mri.dim]));
disp(CC)


% Get centroid(s) of region(s)
centroids = nan(CC.NumObjects, 3);
for icen = 1:CC.NumObjects
    centroids(icen,:) = mean(interp.pos(CC.PixelIdxList{icen},:));
    fprintf('Centroid %i:\nx = %f\ny = %f\nz = %f\n\n', ...
        icen, centroids(icen,1), centroids(icen,2), centroids(icen,3))
end


% For each leadfield voxel inside the brain,
% simply compute the distance to the centroid(s), 
% selecting all voxels within 1cm and inside the brain.

% Very broadly narrow down voxels for which to compute this
lf_mask_right = leadfield.pos(:,1) >=0;
lf_mask_posterior = leadfield.pos(:,2) <= -20;
lf_mask = lf_mask_right & lf_mask_posterior & leadfield.inside;

% Get overview of leadfield voxel distaneces to centroids
figure; 
pos_candidates = leadfield.pos(lf_mask,:);
lf_vox_for_centroid = cell(CC.NumObjects, 1);
for icen = 1:CC.NumObjects
    centroid = centroids(icen,:);
    eucl_dist = nan(sum(lf_mask),1);
    for ivox = 1:sum(lf_mask)
        pos = pos_candidates(ivox,:);
        dist_mat = dist([pos; centroid]');
        eucl_dist(ivox) = dist_mat(1,2);
    end
    subplot(3,4,icen)
    hist(eucl_dist, 100); hold on
    xlabel('euclidean distance (mm)'); ylabel('count')
    title(sprintf('distance from voxels to centroid %i', icen))
    
    % Select leadfield voxels within 10mm of centroid
    min_dist = min(eucl_dist);
    voxoi = pos_candidates(eucl_dist < 10,:);
    disp(voxoi)
    lf_vox_for_centroid{icen} = voxoi;
end


% Visualize ROI for fslide source analysis (low quality)
collim = max(abs(min(dummy.effect)), abs(max(dummy.effect)));
cfg = [];
cfg.location      = centroids;
cfg.method        = 'ortho';
cfg.funparameter  = 'effect';
%cfg.maskparameter = 'mask_effect';
cfg.funcolorlim   = [-collim, collim];
cfg.atlas = vtpm;
ft_sourceplot(cfg, interp);
colormap('viridis')
saveas(gcf, fullfile(fig_save_dir, 'ortho_slices_at_centroid_of_roi_low_qual.png'))


% Visualize ROI for fslide source analysis (high quality)
useopengl = false;
set(0, 'defaultFigureRenderer', 'painters')

collim = max(abs(min(dummy.effect)), abs(max(dummy.effect)));
cfg = [];
cfg.location      = centroids;
cfg.method        = 'ortho';
cfg.funparameter  = 'effect';
%cfg.maskparameter = 'mask_effect';
cfg.funcolorlim   = [-collim, collim];
cfg.atlas = vtpm;
ft_sourceplot(cfg, interp);
colormap('viridis')
saveas(gcf, fullfile(fig_save_dir, 'ortho_slices_at_centroid_of_roi.svg'))


%% Save

save(fullfile(src_dir, 'ROI_1f2f_vs_all_noise_in_coi.mat'), ...
    'lf_vox_for_centroid', 'centroids');


% eof

