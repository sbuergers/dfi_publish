% Figures for time collapsed analysis of dprime and bias
% for high and low alpha frequency quantiles for either high or low
% alpha power trials (select manually which one!). 
%
% Parent script(s): 
%   dfi_fslide_taking_sessions_into_account_prep.m
%   dfi_fslide_taking_sdtparams_hi_v_lo_power.m
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
% Last modified May. 2021


clc
close all
clear all


%% SELECT: hi_pow, lo_pow

pow_cond = 'lo_pow';

% Ran hi_pow on 16.05.2021, 11:31h
% Ran lo_pow on 16.05.2021, 11:30h


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
main_dir = fullfile(data_dir, 'sdt', 'freq_slide', pow_cond);
main_save_dir = fullfile('dfi_experiment_figures', 'Paper_figures', 'iAF', ...
                         'iAF_within', 'sd_params', pow_cond);
mkdir(main_save_dir);

% Colours for plotting
col1 = {[1 0 0], [1 0 0], [1 0 0], [1 0 0]};
col2 = {[0 0 1], [0 0 1], [0 0 1], [0 0 1]};

% use opengl
useopengl = false;

% prepare data
condvect = {'v2', 'Fus', 'Fis'}; 


%% Figure 1 (Instantaneous freq at iAF - dprime)

fh1 = figure('color', [1 1 1], 'Position', [0, 0, 427, 705]);
ha = tight_subplot(6, 4, [0.02 0.02], [0.02], [0.02]);


%% Column 1 (yn_intermsoas, sensor)
clear fslide_GA tif sem_f1 sem_f2 stat

load(fullfile(main_dir, 'sd_params_d_c_yesno_tcollapse.mat'));

axes(ha(1));

xl = [0.5 3.5];
yl = [0.9, 2.8];

for icond = 1:3

    errorbar((icond-1)+0.9, mean(dp_mat(:,icond,1)), dp_within_SE(icond,1), dp_within_SE(icond,1), ...
        'Color', col1{icond}, 'LineWidth', 1.5); hold on
    errorbar((icond-1)+1.1, mean(dp_mat(:,icond,3)), dp_within_SE(icond,3), dp_within_SE(icond,3), ...
        'Color', col2{icond}, 'LineWidth', 1.5);
    plot((icond-1)+0.9, mean(dp_mat(:,icond,1)), ...
        'ko', 'MarkerSize', 3, 'MarkerFaceColor', col1{icond});
    plot((icond-1)+1.1, mean(dp_mat(:,icond,3)), ...
        'ko', 'MarkerSize', 3, 'MarkerFaceColor', col2{icond});
    
    xlim(xl)
    ylim(yl)
    yticks = 0:0.5:100;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    set(gca,'XColor','k','YColor','k')
    
end % condition loop


%% Column 2 (ynt sensor)
clear fslide_GA tif sem_f1 sem_f2 stat

load(fullfile(main_dir, 'sd_params_d_c_tcollapse.mat'));

axes(ha(2));

xl = [0.5 3.5];
yl = [0.9, 2.8];

for icond = 1:3
    
    errorbar((icond-1)+0.9, mean(dp_mat(:,icond,1)), dp_within_SE(icond,1), dp_within_SE(icond,1), ...
        'Color', col1{icond}, 'LineWidth', 1.5); hold on
    errorbar((icond-1)+1.1, mean(dp_mat(:,icond,3)), dp_within_SE(icond,3), dp_within_SE(icond,3), ...
        'Color', col2{icond}, 'LineWidth', 1.5);
    plot((icond-1)+0.9, mean(dp_mat(:,icond,1)), ...
        'ko', 'MarkerSize', 3, 'MarkerFaceColor', col1{icond});
    plot((icond-1)+1.1, mean(dp_mat(:,icond,3)), ...
        'ko', 'MarkerSize', 3, 'MarkerFaceColor', col2{icond});
    
    xlim(xl)
    ylim(yl)
    yticks = 0:0.5:100;  set(gca, 'YTick', yticks);
    set(gca, 'xticklabel', [])
    set(gca, 'yticklabel', [])
    
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
    set(gca,'XColor','k','YColor','k')
    
end % condition loop

fh.Renderer = 'painters'; 
saveas(fh1,fullfile(main_save_dir, 'fslide_dprime_tcollapse_svg.svg'))
close all


%% Figure 2 (Bayes factors - dprime)

fh2 = figure('color', [1 1 1], 'Position', [0, 0, 427, 400]);
ha = tight_subplot(6, 4,[0.02 0.02],[0.02],[0.02]);

fprintf('\n+++++ STATISTICS DPRIME +++++\n')


%% Column 1 (yn_intermsoas, sensor)

load(fullfile(main_dir, 'sd_params_d_c_yesno_tcollapse.mat'));

axes(ha(1));

xl = [0.5 3.5];
yl = [-1, 1];

[bfs, tvals, pvals, nobs] = deal(nan(3,1));
fprintf('\n\n----- yn_intermsoas, sensor ------\n')
for icond = 1:3
    [~, p, ~, stats] = ttest(dp_mat(:,icond,1) - dp_mat(:,icond,3));
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

load(fullfile(main_dir, 'sd_params_d_c_tcollapse.mat'));

axes(ha(2));

xl = [0.5 3.5];
yl = [-1, 1];

[bfs, tvals, pvals, nobs] = deal(nan(3,1));
fprintf('\n\n----- ynt sensor ------\n')
for icond = 1:3
    [~, p, ~, stats] = ttest(dp_mat(:,icond,1) - dp_mat(:,icond,3));
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

% Save
fh2.Renderer = 'painters';
saveas(fh2, fullfile(main_save_dir, 'fslide_dprime_tcollapse_bf_svg.svg'))



%% Figure 3 (Instantaneous freq at iAF - criterion)

fh1 = figure('color', [1 1 1], 'Position', [0, 0, 427, 705]);
ha = tight_subplot(6, 4,[0.02 0.02],[0.02],[0.02]);


%% Column 1 (yn_intermsoas, sensor)

load(fullfile(main_dir, 'sd_params_d_c_yesno_tcollapse.mat'));

axes(ha(1));

xl = [0.5 3.5];
yl = [-1.25, 1.25];

for icond = 1:3

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

load(fullfile(main_dir, 'sd_params_d_c_tcollapse.mat'));

axes(ha(2));

xl = [0.5 3.5];
yl = [-1.25, 1.25];

for icond = 1:3
    
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

% Save
fh.Renderer = 'painters'; 
saveas(fh1,fullfile(main_save_dir, 'fslide_criterion_tcollapse_dfi_bug_fixed.svg'))
close all


%% Figure 4 (Bayes factors - criterion)

fh2 = figure('color', [1 1 1], 'Position', [0, 0, 427, 400]);
ha = tight_subplot(6, 4,[0.02 0.02],[0.02],[0.02]);

fprintf('\n+++++ STATISTICS CRITERION +++++\n')


%% Column 1 (yn_intermsoas, sensor)

axes(ha(1));

load(fullfile(main_dir, 'sd_params_d_c_yesno_tcollapse.mat'));

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

load(fullfile(main_dir, 'sd_params_d_c_tcollapse.mat'));

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

% Save
fh2.Renderer = 'painters'; 
saveas(fh2, fullfile(main_save_dir, 'fslide_criterion_tcollapse_bf_svg.svg'))
close all


% eof

