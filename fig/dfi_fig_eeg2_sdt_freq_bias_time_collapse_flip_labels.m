% Figures for time collapsed analysis of dprime and bias
% for high and low alpha frequency quantiles after flipping labels based
% on other yes-no experiment.
%
% Parent script(s): 
%   dfi_fslide_yesno_taking_sessions_into_account_sdtparams.m
%   dfi_fslide_taking_sessions_into_account_sdtparams.m
%
%
% DETAILS
%
% After testing for correlations of alpha frequency effects wrt either bias
% or sensitivity over tasks (yes-no, yes-no threshold) within subjects, 
% it is prudent to do a follow up analysis where those correlations were
% significant.
%
% A significant correlation in the absence of effects for single
% experiments indiciates that some subjects might have an effect exactly
% opposite to other subjects (this would mean also opposite to the alpha
% temporal resolution hypothesis).
%
% Now, we want to assess whether an effect emerges after flipping the
% labels for high and low alpha frequency terciles for each subject based
% on their effect in the other experiment, if there was a significant 
% correlation of effects over experiments within subjects.
%
% Specifically, effects that had a significant pearson correlation were:
% * Bias_centre, 1 Sound: N=20, r=0.569, p=0.0088, rho=0.206, p=0.382
% * Bias_centre, 2 Sounds: N=20, r=0.529, p=0.0164, rho=0.5233, p=0.0193
% Though it is noteworthy, that the effect for 1 Sound is driven by a
% single observation/outlier, and is far from significant when Spearman's
% Rho is computed.
% 
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


clc
close all
clear all

% experiment script folder
try
    addpath(genpath('dfi'))
catch
    warning('Cannot find dfi folder')
end

% run startup function
dfi_startup

% paths
data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
main_dir = fullfile(data_dir, 'sdt', 'freq_slide');
main_save_dir = fullfile('dfi_experiment_figures', 'Paper_figures', 'iAF', ...
                         'iAF_within', 'sd_params', 'flip_sign');
mkdir(main_save_dir);
                    
% Condition labels
condvect = {'v2', 'Fus', 'Fis'}; 

% Colours for plotting
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]};
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]};

% use opengl
useopengl = false;


%% Prepare data

load(fullfile(main_dir, 'sd_params_d_c_yesno_tcollapse.mat'));
yn_c_mat = c_mat; 
clear c_mat c_within_SE

load(fullfile(main_dir, 'sd_params_d_c_tcollapse.mat'));
ynt_c_mat = c_mat;
clear c_mat c_within_SE

% alpha frequency effect on bias for yes-no and yes-no threshold exp.
eff_yn = squeeze(yn_c_mat(:,:,1)-yn_c_mat(:,:,3));
eff_ynt = squeeze(ynt_c_mat(:,:,1)-ynt_c_mat(:,:,3));

% we want to flip the sign of the yes-no effect based on the sign of the
% yes-no threshold effect for each participant (we can also do it vice
% versa). So get the sign for both yn and ynt:
sign_yn = sign(eff_yn);
sign_ynt = sign(eff_ynt);

% Look at effect of sign flipping
% [eff_yn, sign_ynt, eff_yn .* sign_ynt]


%% Figure 3 (Instantaneous freq at iAF - criterion)

fh1 = figure('color', [1 1 1], 'Position', [0, 0, 427, 705]);
ha = tight_subplot(6, 4,[0.02 0.02],[0.02],[0.02]);


%% Column 1 (yn_intermsoas, sensor)
% --> Flip sign based on yes-no threshold

c_mat = sign_ynt .* yn_c_mat;

%[squeeze(c_mat(:,:,1)-c_mat(:,:,3)), eff_yn .* sign_ynt]

axes(ha(1));

xl = [0.5 3.5];
yl = [-0.75, 0.75];

c_within_SE = nan(3,3);
for icond = 1:3

    data = squeeze(c_mat(:,icond,:));
    subj_mean = nanmean(data,2);
    grand_mean = nanmean(subj_mean,1);
    data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    c_within_SE(icond,:) = nanstd(data_corr,0,1)./sqrt(size(data_corr,1));

    errorbar((icond-1)+0.9, mean(c_mat(:,icond,1)), c_within_SE(icond,1), c_within_SE(icond,1), ...
        'Color', col1{icond}, 'LineWidth', 1.5); hold on
    errorbar((icond-1)+1.1, mean(c_mat(:,icond,3)), c_within_SE(icond,3), c_within_SE(icond,3), ...
        'Color', col2{icond}, 'LineWidth', 1.5);
    plot((icond-1)+0.9, mean(c_mat(:,icond,1)), ...
        'ko', 'MarkerSize', 3, 'MarkerFaceColor', col1{icond});
    plot((icond-1)+1.1, mean(c_mat(:,icond,3)), ...
        'ko', 'MarkerSize', 3, 'MarkerFaceColor', col2{icond});
    
    xlim(xl)
    ylim(yl)
    yticks = -100:0.5:100;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    set(gca,'XColor','k','YColor','k')
    
end % condition loop


%% Column 2 (ynt sensor)
% --> Flip sign based on yes-no

c_mat = sign_yn .* ynt_c_mat;

axes(ha(2));

xl = [0.5 3.5];
yl = [-0.75, 0.75];

c_within_SE = nan(3,3);
for icond = 1:3

    data = squeeze(c_mat(:,icond,:));
    subj_mean = nanmean(data,2);
    grand_mean = nanmean(subj_mean,1);
    data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [1,size(data,2),size(data,3)]);
    c_within_SE(icond,:) = nanstd(data_corr,0,1)./sqrt(size(data_corr,1));
    
    errorbar((icond-1)+0.9, mean(c_mat(:,icond,1)), c_within_SE(icond,1), c_within_SE(icond,1), ...
        'Color', col1{icond}, 'LineWidth', 1.5); hold on
    errorbar((icond-1)+1.1, mean(c_mat(:,icond,3)), c_within_SE(icond,3), c_within_SE(icond,3), ...
        'Color', col2{icond}, 'LineWidth', 1.5);
    plot((icond-1)+0.9, mean(c_mat(:,icond,1)), ...
        'ko', 'MarkerSize', 3, 'MarkerFaceColor', col1{icond});
    plot((icond-1)+1.1, mean(c_mat(:,icond,3)), ...
        'ko', 'MarkerSize', 3, 'MarkerFaceColor', col2{icond});
    
    xlim(xl)
    ylim(yl)
    yticks = -100:0.5:100;;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    set(gca,'XColor','k','YColor','k')
    
end % condition loop

fh.Renderer = 'painters'; 
saveas(fh1,fullfile(main_save_dir, 'fslide_criterion_tcollapse_dfi_bug_fixed.svg'))
close all


%% Figure 4 (Bayes factors - criterion)

fh2 = figure('color', [1 1 1], 'Position', [0, 0, 427, 400]);
ha = tight_subplot(6, 4,[0.02 0.02],[0.02],[0.02]);

fprintf('\n+++++ STATISTICS CRITERION +++++\n')


%% Column 1 (yn_intermsoas, sensor)

axes(ha(1));

c_mat = yn_c_mat;
c_mat(:, :, 1) = sign_ynt .* c_mat(:, :, 1);
c_mat(:, :, 3) = sign_ynt .* c_mat(:, :, 3);

xl = [0.5 3.5];
yl = [-1, 1];

[bfs, tvals, pvals, nobs] = deal(nan(3,1));
fprintf('\n\n----- yn_intermsoas, sensor ------\n')
for icond = 1:3
    [~, p, ~, stats] = ttest(c_mat(:,icond,1) - c_mat(:,icond,3));
    t = stats.tstat;
    df = stats.df;
    bfs(icond) = t1smpbf(t, df+1);
    tvals(icond) = t;
    pvals(icond) = p;
    nobs = df+1;
    
    fprintf('\nCondition %i\n\nt = %f\np = %f\nN = %i\nbf = %f', ...
        icond, t, p, df+1, bfs(icond))
    disp(' ')
end

bar(log10(bfs)); hold on 

xlim(xl)
ylim(yl)
set(gca, 'xticklabel', [])
set(gca, 'yticklabel', [])

box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')


%% Column 2 (ynt sensor)

axes(ha(2));

c_mat = ynt_c_mat;
c_mat(:, :, 1) = sign_yn .* c_mat(:, :, 1);
c_mat(:, :, 3) = sign_yn .* c_mat(:, :, 3);

xl = [0.5 3.5];
yl = [-1, 1];

[bfs, tvals, pvals, nobs] = deal(nan(3,1));
fprintf('\n\n----- ynt sensor ------\n')
for icond = 1:3
    [~, p, ~, stats] = ttest(c_mat(:,icond,1) - c_mat(:,icond,3));
    t = stats.tstat;
    df = stats.df;
    bfs(icond) = t1smpbf(t, df+1);
    tvals(icond) = t;
    pvals(icond) = p;
    nobs = df+1;
    
    fprintf('\nCondition %i\n\nt = %f\np = %f\nN = %i\nbf = %f', ...
        icond, t, p, df+1, bfs(icond))
    disp(' ')
end

bar(log10(bfs)); hold on 

xlim(xl)
ylim(yl)
set(gca, 'xticklabel', [])
set(gca, 'yticklabel', [])

box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')
xl = [-0.6 -0.1];

fh2.Renderer = 'painters'; 
saveas(fh2, fullfile(main_save_dir, 'fslide_criterion_tcollapse_bf_svg.svg'))
close all


% eof

