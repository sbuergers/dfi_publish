% Make correlation plots to assess consistency of PF threshold
% between yes-no and 2IFC task for each sound condition. In addition,
% compare each task's threshold estimate with the staircase estimate
% of the yes-no thresh. task (3 (task comparison) x 3 (# beeps)
% comparisons)
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


% Yes-no threshold
% Import yes-no threshold data and get individual SOAs
load(fullfile('dfi_experiment_data', 'data', 'experiment', 'd701to727_ynt.mat'))
clear d7*

subjects = unique(dall.partid);
trialtypes = [3, 6, 8, 9];
soamat = nan(20,4);
for isubj = 1:20
    for itrl = 1:length(trialtypes)
        soamat(isubj,itrl) = unique(dall.soa(dall.partid == subjects(isubj) & dall.trlid == trialtypes(itrl)) );
    end
end
ynt_staircase_estimate = soamat;


%% Parameter consistency over tasks

% Some settings
col_vect = [[0 0.6 0]; [0 0 1]; [1 0 0]; [0 0 0]];
col_lines = [[0 0.8 0]; [0 0 1]; [1 0 0]; [0 0 0]];
opacity = 0.5; % transparency (alpha)
lw = 0.6;
ls = '-';

xl = [-0.6 0.4]; 
yl = xl;

% Initialize 3 (task comparison) x 3 (condition) Bayes factor matrix
BFs = nan(3,3);
    

%% Scatter plots 1: yesno versus 2ifc
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

for i = 1:3
    BFs(1, i) = corrbf(r(i), N(i));
end

fprintf('\n\n---- Scatter plots 1: yesno versus 2ifc ----')
for icond = 1:3
    fprintf('\n\nCondition %i\n', icond)
    fprintf('\nN = %i', N(icond))
    fprintf('\nr = %f', r(icond))
    fprintf('\nrho = %f', spearRho(icond))
    fprintf('\np = %f', pval(icond))
    fprintf('\nbf = %f', BFs(1, icond))
end
fprintf('\n\n')


%% Scatter plots 2: 2ifc versus ynt staircase estimate

x = ifc_beta_threshold;
y = ynt_staircase_estimate;
    
for ic = 1:3
    
    axes(ha(4+ic));

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

for i = 1:3
    BFs(2, i) = corrbf(r(i), N(i));
end

fprintf('\n\n---- Scatter plots 1: 2ifc vs yesno threshold ----')
for icond = 1:3
    fprintf('\n\nCondition %i\n', icond)
    fprintf('\nN = %i', N(icond))
    fprintf('\nr = %f', r(icond))
    fprintf('\nrho = %f', spearRho(icond))
    fprintf('\np = %f', pval(icond))
    fprintf('\nbf = %f', BFs(2, icond))
end
fprintf('\n\n')


%% Scatter plots 3: yesno versus ynt staircase estimate

x = yn_beta_threshold;
y = ynt_staircase_estimate;
    
for ic = 1:3
    
    axes(ha(8+ic));

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

for i = 1:3
    BFs(3, i) = corrbf(r(i), N(i));
end

fprintf('\n\n---- Scatter plots 3: yesno versus yesno threshold ----')
for icond = 1:3
    fprintf('\n\nCondition %i\n', icond)
    fprintf('\nN = %i', N(icond))
    fprintf('\nr = %f', r(icond))
    fprintf('\nrho = %f', spearRho(icond))
    fprintf('\np = %f', pval(icond))
    fprintf('\nbf = %f', BFs(3, icond))
end
fprintf('\n\n')

% save figure
fh0.Renderer = 'painters'; 
saveas(fh0, fullfile(figdir, 'Consistency_pf_threshold_svg.svg'))
close all


%% Plot bayes factors

fh1 = figure('color', [1 1 1], 'Position', [0, 0, 427, 350]);
ha = tight_subplot(3, 4,[0.01 0.03],[0.02],[0.02]);

xl = [0, 1];
yl = [-1, +1];

for itask = 1:3
    axes(ha(itask));

    line([xl], [0 0])
    line([xl], [0.5, 0.5])
    line([xl], [-0.5, -0.5])
    line([0.25, 0.25], [0, log10(BFs(itask, 1))])
    line([0.5, 0.5], [0, log10(BFs(itask, 2))])
    line([0.75, 0.75], [0, log10(BFs(itask, 3))])
    box off
end

% save figure
fh1.Renderer = 'painters'; 
saveas(fh1, fullfile(figdir, 'Consistency_pf_threshold_BFs_svg.svg'))
close all

 
% eof

