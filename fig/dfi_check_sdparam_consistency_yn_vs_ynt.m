% Check if participants have consistently higher or lower sensitivity or
% bias in both experiments (yesno and yes-nothreshold); irrespective of
% alpha frequency. 
%
% Parent script(s): 
%   dfi_fslide_yesno_taking_sessions_into_account_sdtparams.m
%   dfi_fslide_taking_sessions_into_account_sdtparams.m
%
% Sibling script(s):
%   dfi_check_fslide_effect_consistency_yn_vs_ynt_no_reg.m
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
data_dir = fullfile('dfi_experiment_data', 'data', 'experiment');
main_save_dir = fullfile('dfi_experiment_figures', 'Paper_figures', 'beh');
mkdir(main_save_dir)

% prepare data
condvect = {'v2', 'Fus', 'Fis'}; 

% Load yes-no sd params for frequency
load(fullfile( data_dir, 'SD_params_yn_and_ynt.mat' ));
dp_ynt = dP_ynt';
c_ynt  = C_ynt';
dp_yn = dP_yesno';
c_yn  = C_yesno';


% Some settings
col_vect = [[0 0.6 0]; [0 0 1]; [1 0 0]; [0 0 0]];
col_lines = [[0 0.8 0]; [0 0 1]; [1 0 0]; [0 0 0]];
opacity = 0.5; % transparency (alpha)
lw = 0.6;
ls = '-';

xl = [0 5]; 
yl = xl;


% Scatter plot: Collapse over time, check yesno versus ynt
fh0 = figure('color', [1 1 1], 'Position', [0, 0, 427, 350]);
ha = tight_subplot(3, 4,[0.01 0.03],[0.02],[0.02]);


%% Row 1: Dprime

for icond = 1:3
    ic = icond;
    x = dp_yn;
    y = dp_ynt;
    
    axes(ha(icond));

    xrange = xl(2) - xl(1);
    yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(x(:,ic), y(:,ic), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        N(1,ic) = sum(~isnan(x(:,ic)) & ~isnan(y(:,ic)));
        [spearRho(1,ic), pval(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Spearman', 'rows', 'complete');
        [rx(1,ic), pval_r(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Pearson', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(x(:,ic), y(:,ic), 'one');
        line([min(min(x))  max(max(x))], [b0(1,ic)+b1(1,ic)*min(min(x)) b0(1,ic)+b1(1,ic)*max(max(x))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = -0:1:5; set(gca, 'XTick', xticks)
        yticks = -0:1:5;  set(gca, 'YTick', yticks);
    end
    line([xl], [0 0], 'color', [0.5 0.5 0.5])
    line([0 0], [yl], 'color', [0.5 0.5 0.5])
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
end

yl
N
r
pval_r
spearRho
pval
r_dp = r;
n_dp = N;


%% Row 2: Criterion

xl = [-3 3];
yl = xl;

for icond = 1:3
    ic = icond;
    x = c_yn;
    y = c_ynt;
    
    axes(ha(icond+4));

    xrange = xl(2) - xl(1);
    yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(x(:,ic), y(:,ic), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        N(1,ic) = sum(~isnan(x(:,ic)) & ~isnan(y(:,ic)));
        [spearRho(1,ic), pval(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Spearman', 'rows', 'complete');
        [rx(1,ic), pval_r(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Pearson', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(x(:,ic), y(:,ic), 'one');
        line([min(min(x))  max(max(x))], [b0(1,ic)+b1(1,ic)*min(min(x)) b0(1,ic)+b1(1,ic)*max(max(x))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = -3:1:3; set(gca, 'XTick', xticks)
        yticks = -3:1:3;  set(gca, 'YTick', yticks);
    end
    line([xl], [0 0], 'color', [0.5 0.5 0.5])
    line([0 0], [yl], 'color', [0.5 0.5 0.5])
    box off
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
end

yl
N
r
pval_r
spearRho
pval
r_c = r;
n_c = N;


fh.Renderer = 'painters'; 
saveas(fh0, fullfile(main_save_dir, 'Consistency_sdparams_svg.svg'))
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
saveas(fh1, fullfile(main_save_dir, 'Consistency_sdparams_svg_BFs.svg'))
close all


% eof

