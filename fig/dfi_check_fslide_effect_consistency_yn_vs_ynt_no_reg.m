% Check if participants have consistently higher or lower sensitivity or
% bias depending on alpha frequency in both experiments (yesno and
% yes-nothreshold). 
% Supplementary analysis: No regressing out of nuisance variables.
%
% Parent script(s): 
%   dfi_fslide_yesno_sdtparams.m
%   dfi_fslide_sdtparams.m
%
% Sibling script(s):
%   dfi_check_fslide_effect_consistency_yn_vs_ynt.m
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
% Last modified Mar. 2021


clc
close all
clear all

% paths
data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
main_dir = fullfile(data_dir, 'sdt', 'freq_slide', 'no_regress');
main_save_dir = fullfile('dfi_experiment_figures', 'Paper_figures', 'iAF', ...
                         'iAF_within', 'sd_params', 'no_regress');
src_dir = fullfile(data_dir, 'source_analysis', 'sdt', 'freq_slide', 'no_regress');
src_save_dir = fullfile('dfi_experiment_figures', 'Paper_figures', 'iAF', ...
                        'iAF_within', 'sd_params', 'no_regress');
                    
mkdir(main_save_dir)
mkdir(src_save_dir)

% prepare data
condvect = {'v2', 'Fus', 'Fis'}; 

% Load yes-no sd params for frequency
load(fullfile(main_dir, 'sd_params_d_c_tcollapse.mat'));

dp_ynt = dp_mat; clear dp_mat
c_ynt  = c_mat; clear c_mat

% Load yes-no threshold sd params for frequency
load(fullfile(main_dir, 'sd_params_d_c_yesno_tcollapse.mat'));

dp_yn = dp_mat; clear dp_mat
c_yn  = c_mat; clear c_mat


% Some settings
col_vect = [[0 0.6 0]; [0 0 1]; [1 0 0]; [0 0 0]];
col_lines = [[0 0.8 0]; [0 0 1]; [1 0 0]; [0 0 0]];
opacity = 0.5; % transparency (alpha)
lw = 0.6;
ls = '-';

xl = [-0.6 0.4]; 
yl = xl;


% Scatter plot: Collapse over time, check yesno versus ynt
fh0 = figure('color', [1 1 1], 'Position', [0, 0, 427, 350]);
ha = tight_subplot(3, 4,[0.01 0.03],[0.02],[0.02]);


%% Row 1: Dprime

for icond = 1:3
    ic = icond;
    x = squeeze(dp_yn(:,:,1)-dp_yn(:,:,3));
    y = squeeze(dp_ynt(:,:,1)-dp_ynt(:,:,3));
    
    axes(ha(icond));

    % Outliers are beyond 3 scaled median absolute deviations (MAD)
    % MAD = K * median(|Ai - median(A)|), where
    % i = 1, 2, ..., N; and K ~= 1.4826 (see help isoutlier)
    otl = isoutlier(x) | isoutlier(y);

    xl = [min(min(x(:,ic)))*1.15, max(max(x(:,ic)))*1.15]; xrange = xl(2) - xl(1);
    xlall(:,ic) = xl;
    yl = [min(min(y))*1.15, max(max(y))*1.15]; yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(x(~otl(:, ic),ic), y(~otl(:, ic),ic), ...
        col_vect(icond,:), opacity, dtsz, 25); hold on
    if any(otl(:, ic))
        plot(x(otl(:, ic), ic), y(otl(:, ic), ic), '*k', 'markersize', 3)
    end
    N(1,ic) = sum(~isnan(x(~otl(:, ic),ic)) & ~isnan(y(~otl(:, ic),ic)));
    [spearRho(1,ic), pval(1,ic)] = corr(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'type', 'Spearman', 'rows', 'complete');
    [rx(1,ic), pval_r(1,ic)] = corr(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'type', 'Pearson', 'rows', 'complete');
    [r(1,ic),b1(1,ic),b0(1,ic)] = regression(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'one');
    line([min(min(x))  max(max(x))], [b0(1,ic)+b1(1,ic)*min(min(x)) b0(1,ic)+b1(1,ic)*max(max(x))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
    xlim(xl); ylim(yl);
    xticks = -2:0.4:2; set(gca, 'XTick', xticks)
    yticks = -2:0.4:2;  set(gca, 'YTick', yticks);
    line([xl], [0 0], 'color', [0.5 0.5 0.5])
    line([0 0], [yl], 'color', [0.5 0.5 0.5])
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
end

xlall'
yl
N
r
pval_r
spearRho
pval
r_dp = r;
n_dp = N;


%% Row 2: Criterion

for icond = 1:3
    ic = icond;
    x = squeeze(c_yn(:,:,1)-c_yn(:,:,3));
    y = squeeze(c_ynt(:,:,1)-c_ynt(:,:,3));
    
    axes(ha(icond+4));

    % Outliers are beyond 3 scaled median absolute deviations (MAD)
    % MAD = K * median(|Ai - median(A)|), where
    % i = 1, 2, ..., N; and K ~= 1.4826 (see help isoutlier)
    otl = isoutlier(x) | isoutlier(y);

    xl = [min(min(x(:,ic)))*1.15, max(max(x(:,ic)))*1.15]; xrange = xl(2) - xl(1);
    xlall(:,ic) = xl;
    yl = [min(min(y))*1.15, max(max(y))*1.15]; yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(x(~otl(:, ic),ic), y(~otl(:, ic),ic), ...
        col_vect(icond,:), opacity, dtsz, 25); hold on
    if any(otl(:, ic))
        plot(x(otl(:, ic), ic), y(otl(:, ic), ic), '*k', 'markersize', 3)
    end
    N(1,ic) = sum(~isnan(x(~otl(:, ic),ic)) & ~isnan(y(~otl(:, ic),ic)));
    [spearRho(1,ic), pval(1,ic)] = corr(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'type', 'Spearman', 'rows', 'complete');
    [rx(1,ic), pval_r(1,ic)] = corr(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'type', 'Pearson', 'rows', 'complete');
    [r(1,ic),b1(1,ic),b0(1,ic)] = regression(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'one');
    line([min(min(x))  max(max(x))], [b0(1,ic)+b1(1,ic)*min(min(x)) b0(1,ic)+b1(1,ic)*max(max(x))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
    xlim(xl); ylim(yl);
    xticks = -2:0.2:2; set(gca, 'XTick', xticks)
    yticks = -2:0.2:2;  set(gca, 'YTick', yticks);
    line([xl], [0 0], 'color', [0.5 0.5 0.5])
    line([0 0], [yl], 'color', [0.5 0.5 0.5])
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
end

xlall'
yl
N
r
pval_r
spearRho
pval
r_c = r;
n_c = N;


fh.Renderer = 'painters'; 
saveas(fh0, fullfile(main_save_dir, 'Consistency_freq_sdparams_svg.svg'))
close all


%% Bayes factor figures

[BFs_dp, BFs_c] = deal(nan(3,1));
for i = 1:3
    BFs_dp(i) = corrbf(r_dp(i), n_dp(i));
    BFs_c(i) = corrbf(r_c(i), n_c(i));
end

fh1 = figure('color', [1 1 1], 'Position', [0, 0, 427, 350]);
ha = tight_subplot(3, 4,[0.01 0.03],[0.02],[0.02]);

%% 1: Dprime
    
xl = [0, 1];
yl = [-1, +1];

axes(ha(1));

line([xl], [0 0])
line([xl], [0.5, 0.5])
line([xl], [-0.5, -0.5])
line([0.25, 0.25], [0, log10(BFs_dp(1))])
line([0.5, 0.5], [0, log10(BFs_dp(2))])
line([0.75, 0.75], [0, log10(BFs_dp(3))])
box off


%% 2: Bias

axes(ha(2))

line([xl], [0 0])
line([xl], [0.5, 0.5])
line([xl], [-0.5, -0.5])
line([0.25, 0.25], [0, log10(BFs_c(1))])
line([0.5, 0.5], [0, log10(BFs_c(2))])
line([0.75, 0.75], [0, log10(BFs_c(3))])
box off

fh.Renderer = 'painters'; 
saveas(fh1, fullfile(main_save_dir, 'Consistency_freq_sdparams_svg_BFs.svg'))
close all



%% LCMV (same analysis in source space)

% Load yes-no sd params for frequency
load(fullfile(src_dir, 'sd_params_d_c_tcollapse.mat'));

dp_ynt = dp_mat; clear dp_mat
c_ynt  = c_mat; clear c_mat

% Load yes-no threshold sd params for frequency
load(fullfile(src_dir, 'sd_params_d_c_yesno_tcollapse.mat'));

dp_yn = dp_mat; clear dp_mat
c_yn  = c_mat; clear c_mat


% Some settings
col_vect = [[0 0.6 0]; [0 0 1]; [1 0 0]; [0 0 0]];
col_vect_ci = [[0.75 1 0.75]; [0.75 0.75 1]; [1 0.75 0.75]; [0.75 0.75 0.75]];
col_lines = [[0 0.8 0]; [0 0 1]; [1 0 0]; [0 0 0]];
opacity = 0.5; % transparency (alpha)
lw = 0.6;
ls = '-';

xl = [-0.6 0.4]; 
yl = xl;
plot_ci = false;


% Scatter plot: Collapse over time, check yesno versus ynt
fh0 = figure('color', [1 1 1], 'Position', [0, 0, 427, 350]);
ha = tight_subplot(3, 4,[0.01 0.03],[0.02],[0.02]);


%% Row 1: Dprime
for icond = 1:3
    ic = icond;
    x = squeeze(dp_yn(:,:,1)-dp_yn(:,:,3));
    y = squeeze(dp_ynt(:,:,1)-dp_ynt(:,:,3));
    
    axes(ha(icond));

    % Outliers are beyond 3 scaled median absolute deviations (MAD)
    % MAD = K * median(|Ai - median(A)|), where
    % i = 1, 2, ..., N; and K ~= 1.4826 (see help isoutlier)
    otl = isoutlier(x) | isoutlier(y);

    xl = [min(min(x(:,ic)))*1.15, max(max(x(:,ic)))*1.15]; xrange = xl(2) - xl(1);
    xlall(:,ic) = xl;
    yl = [min(min(y))*1.15, max(max(y))*1.15]; yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(x(~otl(:, ic),ic), y(~otl(:, ic),ic), ...
        col_vect(icond,:), opacity, dtsz, 25); hold on
    if any(otl(:, ic))
        plot(x(otl(:, ic), ic), y(otl(:, ic), ic), '*k', 'markersize', 3)
    end
    N(1,ic) = sum(~isnan(x(~otl(:, ic),ic)) & ~isnan(y(~otl(:, ic),ic)));
    [spearRho(1,ic), pval(1,ic)] = corr(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'type', 'Spearman', 'rows', 'complete');
    [rx(1,ic), pval_r(1,ic)] = corr(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'type', 'Pearson', 'rows', 'complete');
    [r(1,ic),b1(1,ic),b0(1,ic)] = regression(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'one');
    line([min(min(x))  max(max(x))], [b0(1,ic)+b1(1,ic)*min(min(x)) b0(1,ic)+b1(1,ic)*max(max(x))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
    xlim(xl); ylim(yl);
    xticks = -2:0.4:2; set(gca, 'XTick', xticks)
    yticks = -2:0.4:2;  set(gca, 'YTick', yticks);
    line([xl], [0 0], 'color', [0.5 0.5 0.5])
    line([0 0], [yl], 'color', [0.5 0.5 0.5])
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
end

xlall'
yl
N
r
pval_r
spearRho
pval
r_dp = r;
n_dp = N;


%% Row 2: Criterion
for icond = 1:3
    ic = icond;
    x = squeeze(c_yn(:,:,1)-c_yn(:,:,3));
    y = squeeze(c_ynt(:,:,1)-c_ynt(:,:,3));
    
    axes(ha(icond+4));

    % Outliers are beyond 3 scaled median absolute deviations (MAD)
    % MAD = K * median(|Ai - median(A)|), where
    % i = 1, 2, ..., N; and K ~= 1.4826 (see help isoutlier)
    otl = isoutlier(x) | isoutlier(y);

    xl = [min(min(x(:,ic)))*1.15, max(max(x(:,ic)))*1.15]; xrange = xl(2) - xl(1);
    xlall(:,ic) = xl;
    yl = [min(min(y))*1.15, max(max(y))*1.15]; yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(x(~otl(:, ic),ic), y(~otl(:, ic),ic), ...
        col_vect(icond,:), opacity, dtsz, 25); hold on
    if any(otl(:, ic))
        plot(x(otl(:, ic), ic), y(otl(:, ic), ic), '*k', 'markersize', 3)
    end
    N(1,ic) = sum(~isnan(x(~otl(:, ic),ic)) & ~isnan(y(~otl(:, ic),ic)));
    [spearRho(1,ic), pval(1,ic)] = corr(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'type', 'Spearman', 'rows', 'complete');
    [rx(1,ic), pval_r(1,ic)] = corr(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'type', 'Pearson', 'rows', 'complete');
    [r(1,ic),b1(1,ic),b0(1,ic)] = regression(x(~otl(:, ic),ic), y(~otl(:, ic),ic), 'one');
    line([min(min(x))  max(max(x))], [b0(1,ic)+b1(1,ic)*min(min(x)) b0(1,ic)+b1(1,ic)*max(max(x))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
    xlim(xl); ylim(yl);
    xticks = -2:0.2:2; set(gca, 'XTick', xticks)
    yticks = -2:0.2:2;  set(gca, 'YTick', yticks);
    line([xl], [0 0], 'color', [0.5 0.5 0.5])
    line([0 0], [yl], 'color', [0.5 0.5 0.5])
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
end

xlall'
yl
N
r
pval_r
spearRho
pval
r_c = r;
r_n = N;


fh.Renderer = 'painters'; 
saveas(fh0, fullfile(src_save_dir, 'Consistency_freq_sdparams_lcmv_svg.svg'))
close all



%% Bayes factor figures

[BFs_dp, BFs_c] = deal(nan(3,1));
for i = 1:3
    BFs_dp(i) = corrbf(r_dp(i), n_dp(i));
    BFs_c(i) = corrbf(r_c(i), n_c(i));
end

fh1 = figure('color', [1 1 1], 'Position', [0, 0, 427, 350]);
ha = tight_subplot(3, 4,[0.01 0.03],[0.02],[0.02]);


%% 1: Dprime
    
xl = [0, 1];
yl = [-1, +1];

axes(ha(1));

line([xl], [0 0])
line([xl], [0.5, 0.5])
line([xl], [-0.5, -0.5])
line([0.25, 0.25], [0, log10(BFs_dp(1))])
line([0.5, 0.5], [0, log10(BFs_dp(2))])
line([0.75, 0.75], [0, log10(BFs_dp(3))])
box off


%% 2: Bias

axes(ha(2))

line([xl], [0 0])
line([xl], [0.5, 0.5])
line([xl], [-0.5, -0.5])
line([0.25, 0.25], [0, log10(BFs_c(1))])
line([0.5, 0.5], [0, log10(BFs_c(2))])
line([0.75, 0.75], [0, log10(BFs_c(3))])
box off

fh.Renderer = 'painters'; 
saveas(fh1, fullfile(src_save_dir, 'Consistency_freq_sdparams_svg_lcmv_BFs.svg'))
close all


% eof

