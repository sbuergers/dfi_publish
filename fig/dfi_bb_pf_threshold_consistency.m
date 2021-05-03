% Make correlation plots to assess consistency of PF threshold
% between yes-no and 2IFC task for each sound condition.
%
% Parent script(s): 
%   dfi_beta_binom_combine_joined_and_sep_fits.m
%   
% Children script(s): 
%   None
%
% Sibling script(s):
%   dfi_bb_pf_param_consistency_over_tasks_and_conds.m
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


% add project directory to path
addpath(genpath('dfi'))

figdir = fullfile('dfi_experiment_figures', 'PFs', 'beta_binom_weibull');


%% Load data

% 2IFC
% Import beta binomial parameter estimates
folder = '2ifc';
load(fullfile(figdir, folder, 'dataPF_joined_and_sep_fits_combined.mat'), ...
    'threshold_matrix', 'slope_matrix', 'guess_matrix', 'lapse_matrix', 'eta_matrix')
close all

% Also save data in same format as used previously, so other scripts are
% easy to adapt!
pffit_2ifc = cell(1,20);
for isubj = 1:20
    for icond = 1:3
        pffit_2ifc{isubj}{icond}.cond = icond;
        pffit_2ifc{isubj}{icond}.label = 'alpha, beta, gamma, lambda';
        pffit_2ifc{isubj}{icond}.par = [threshold_matrix(isubj,icond), ...
                                        slope_matrix(isubj,icond), ...
                                        guess_matrix(isubj,icond), ...
                                        lapse_matrix(isubj,icond)];                                   
    end
end

ifc_beta_threshold = threshold_matrix;
ifc_beta_slope = slope_matrix;
ifc_beta_guess = guess_matrix;
ifc_beta_lapse = lapse_matrix;


% Yes-no
% Import beta binomial parameter estimates
folder = 'yn_pooled';
load(fullfile(figdir, folder, 'dataPF_joined_and_sep_fits_combined.mat'), ...
    'threshold_matrix', 'slope_matrix', 'guess_matrix', 'lapse_matrix', 'eta_matrix')
close all

% Also save data in same format as used previously, so other scripts are
% easy to adapt!
pffit_yn = cell(1,20);
for isubj = 1:20
    for icond = 1:3
        pffit_yn{isubj}{icond}.cond = icond;
        pffit_yn{isubj}{icond}.label = 'alpha, beta, gamma, lambda';
        pffit_yn{isubj}{icond}.par = [threshold_matrix(isubj,icond), ...
                                      slope_matrix(isubj,icond), ...
                                      guess_matrix(isubj,icond), ...
                                      lapse_matrix(isubj,icond)];                                   
    end
end

yn_beta_threshold = threshold_matrix;
yn_beta_slope = slope_matrix;
yn_beta_guess = guess_matrix;
yn_beta_lapse = lapse_matrix;
yn_beta_eta = eta_matrix;


%% Parameter consistency over tasks

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
ha = tight_subplot(3, 4, [0.01 0.03], [0.02], [0.02]);

x = ifc_beta_threshold;
y = yn_beta_threshold;
    
for ic = 1:3
    
    axes(ha(ic));

    xl = [min(min(x(:,ic)))*0.85, max(max(x(:,ic)))*1.15]; xrange = xl(2) - xl(1);
    xlall(:,ic) = xl;
    yl = [min(min(y(:,ic)))*0.85, max(max(y(:,ic)))*1.15]; yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(x(:,ic), y(:,ic), col_vect(ic,:), opacity, dtsz, 25); hold on
    N(1,ic) = sum(~isnan(x(:,ic)) & ~isnan(y(:,ic)));
    [spearRho(1,ic), pval(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Spearman', 'rows', 'complete');
    [rx(1,ic), pval_r(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Pearson', 'rows', 'complete');
    [r(1,ic),b1(1,ic),b0(1,ic)] = regression(x(:,ic), y(:,ic), 'one');
    line([min(min(x))  max(max(x))], [b0(1,ic)+b1(1,ic)*min(min(x)) b0(1,ic)+b1(1,ic)*max(max(x))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
    xlim(xl); ylim(yl);
    xticks = [0.025, 0.042, 0.05, 0.058, 0.075, 0.108, 0.158, 0.225]; set(gca, 'XTick', xticks)
    yticks = xticks;  set(gca, 'YTick', yticks);
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

% save figure
fh0.Renderer = 'painters'; 
saveas(fh0, fullfile(figdir, 'Consistency_pf_threshold_svg.svg'))
close all


%% Bayes factor figures

BFs = nan(3,1);
for i = 1:3
    BFs(i) = corrbf(r(i), N(i));
end

fh1 = figure('color', [1 1 1], 'Position', [0, 0, 427, 350]);
ha = tight_subplot(3, 4,[0.01 0.03],[0.02],[0.02]);

xl = [0, 1];
yl = [-1, +1];

axes(ha(1));

line([xl], [0 0])
line([xl], [0.5, 0.5])
line([xl], [-0.5, -0.5])
line([0.25, 0.25], [0, log10(BFs(1))])
line([0.5, 0.5], [0, log10(BFs(2))])
line([0.75, 0.75], [0, log10(BFs(3))])
box off

% save figure
fh1.Renderer = 'painters'; 
saveas(fh1, fullfile(figdir, 'Consistency_pf_threshold_BFs_svg.svg'))
close all


%% TODO: Add correlations between 2IFC and Yes-no and staircase

           
% eof

