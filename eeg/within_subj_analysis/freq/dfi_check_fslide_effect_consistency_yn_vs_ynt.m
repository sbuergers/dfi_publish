% Check if participants have consistently higher or lower sensitivity or
% bias depending on alpha frequency in both experiments (yesno and
% yes-nothreshold). 


clc
close all
clear all

% prepare data
condvect = {'v2', 'Fus', 'Fis'}; 

% Load yes-no sd params for frequency
load(fullfile('dfi_experiment_data', 'eeg_data', 'experiment', 'sdt', 'freq_slide', ...
    'sd_params_d_c.mat'));

dp_ynt = dp_mat_cont; clear dp_mat_cont
c_ynt  = c_mat_cont; clear c_mat_cont

% Load yes-no threshold sd params for frequency
load(fullfile('dfi_experiment_data', 'eeg_data', 'experiment', 'sdt', 'freq_slide', ...
    'sd_params_d_c_yesno.mat'));

% collapse over intermediate 4 SOAs (0.05 to 0.108)
dp_mat_cont_full = dp_mat_cont; clear dp_mat_cont
dp_mat_cont = squeeze(nanmean(dp_mat_cont_full(3:6,:,:,:,:)));

c_mat_cont_full = c_mat_cont; clear c_mat_cont
c_mat_cont = squeeze(nanmean(c_mat_cont_full(3:6,:,:,:,:)));

dp_yn = dp_mat_cont; clear dp_mat_cont
c_yn  = c_mat_cont; clear c_mat_cont


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
    x = squeeze(nanmean(dp_yn(:,:,:,1)-dp_yn(:,:,:,3),1));
    y = squeeze(nanmean(dp_ynt(:,:,:,1)-dp_ynt(:,:,:,3),1));
    
    axes(ha(icond));

    xl = [min(min(x(:,ic)))*1.15, max(max(x(:,ic)))*1.15]; xrange = xl(2) - xl(1);
    xlall(:,ic) = xl;
    yl = [min(min(y))*1.15, max(max(y))*1.15]; yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(x(:,ic), y(:,ic), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        N(1,ic) = sum(~isnan(x(:,ic)) & ~isnan(y(:,ic)));
        [spearRho(1,ic), pval(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Spearman', 'rows', 'complete');
        [rx(1,ic), pval_r(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Pearson', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(x(:,ic), y(:,ic), 'one');
        line([min(min(x))  max(max(x))], [b0(1,ic)+b1(1,ic)*min(min(x)) b0(1,ic)+b1(1,ic)*max(max(x))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = -2:0.4:2; set(gca, 'XTick', xticks)
        yticks = -2:0.4:2;  set(gca, 'YTick', yticks);
    end
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


%% Row 2: Criterion

for icond = 1:3
    ic = icond;
    x = squeeze(nanmean(c_yn(:,:,:,1)-c_yn(:,:,:,3),1));
    y = squeeze(nanmean(c_ynt(:,:,:,1)-c_ynt(:,:,:,3),1));
    
    axes(ha(icond+4));

    xl = [min(min(x(:,ic)))*1.15, max(max(x(:,ic)))*1.15]; xrange = xl(2) - xl(1);
    xlall(:,ic) = xl;
    yl = [min(min(y))*1.15, max(max(y))*1.15]; yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(x(:,ic), y(:,ic), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        N(1,ic) = sum(~isnan(x(:,ic)) & ~isnan(y(:,ic)));
        [spearRho(1,ic), pval(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Spearman', 'rows', 'complete');
        [rx(1,ic), pval_r(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Pearson', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(x(:,ic), y(:,ic), 'one');
        line([min(min(x))  max(max(x))], [b0(1,ic)+b1(1,ic)*min(min(x)) b0(1,ic)+b1(1,ic)*max(max(x))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = -2:0.2:2; set(gca, 'XTick', xticks)
        yticks = -2:0.2:2;  set(gca, 'YTick', yticks);
    end
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


fh.Renderer = 'painters'; 
saveas(fh0, fullfile('dfi_experiment_figures', 'Paper_figures', 'iAF', ...
    'iAF_within', 'sd_params', 'Consistency_freq_sdparams_svg.svg'))
close all


%% Bayes factor figures

% Fill in by hand!!!
BFs_dp = [0.3432, 0.1707, 0.7716];
BFs_c = [0.3722, 9.2528, 0.5945];

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
saveas(fh1, fullfile('dfi_experiment_figures', 'Paper_figures', 'iAF', ...
    'iAF_within', 'sd_params', 'Consistency_freq_sdparams_svg_BFs.svg'))
close all


%% LCMV (same analysis in source space)

% Load yes-no sd params for frequency
load(fullfile('dfi_experiment_data', 'eeg_data', 'experiment', 'source_analysis', 'sdt', 'freq_slide', ...
    'sd_params_d_c.mat'));

dp_ynt = dp_mat_cont; clear dp_mat_cont
c_ynt  = c_mat_cont; clear c_mat_cont

% Load yes-no threshold sd params for frequency
load(fullfile('dfi_experiment_data', 'eeg_data', 'experiment', 'source_analysis', 'sdt', 'freq_slide', ...
    'sd_params_d_c_yesno.mat'));

% collapse over intermediate 4 SOAs (0.05 to 0.108)
dp_mat_cont_full = dp_mat_cont; clear dp_mat_cont
dp_mat_cont = squeeze(nanmean(dp_mat_cont_full(3:6,:,:,:,:)));

c_mat_cont_full = c_mat_cont; clear c_mat_cont
c_mat_cont = squeeze(nanmean(c_mat_cont_full(3:6,:,:,:,:)));

dp_yn = dp_mat_cont; clear dp_mat_cont
c_yn  = c_mat_cont; clear c_mat_cont


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
    x = squeeze(nanmean(dp_yn(:,:,:,1)-dp_yn(:,:,:,3),1));
    y = squeeze(nanmean(dp_ynt(:,:,:,1)-dp_ynt(:,:,:,3),1));
    
    axes(ha(icond));

    xl = [min(min(x(:,ic)))*1.15, max(max(x(:,ic)))*1.15]; xrange = xl(2) - xl(1);
    xlall(:,ic) = xl;
    yl = [min(min(y))*1.15, max(max(y))*1.15]; yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(x(:,ic), y(:,ic), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        N(1,ic) = sum(~isnan(x(:,ic)) & ~isnan(y(:,ic)));
        [spearRho(1,ic), pval(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Spearman', 'rows', 'complete');
        [rx(1,ic), pval_r(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Pearson', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(x(:,ic), y(:,ic), 'one');
        line([min(min(x))  max(max(x))], [b0(1,ic)+b1(1,ic)*min(min(x)) b0(1,ic)+b1(1,ic)*max(max(x))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = -2:0.4:2; set(gca, 'XTick', xticks)
        yticks = -2:0.4:2;  set(gca, 'YTick', yticks);
    end
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


%% Row 2: Criterion
for icond = 1:3
    ic = icond;
    x = squeeze(nanmean(c_yn(:,:,:,1)-c_yn(:,:,:,3),1));
    y = squeeze(nanmean(c_ynt(:,:,:,1)-c_ynt(:,:,:,3),1));
    
    axes(ha(icond+4));

    xl = [min(min(x(:,ic)))*1.15, max(max(x(:,ic)))*1.15]; xrange = xl(2) - xl(1);
    xlall(:,ic) = xl;
    yl = [min(min(y))*1.15, max(max(y))*1.15]; yrange = yl(2) - yl(1);
    dtsz = [xrange*.0355,yrange*.0355]; 
    transparentScatter(x(:,ic), y(:,ic), col_vect(icond,:), opacity, dtsz, 25); hold on
    for ic = icond
        N(1,ic) = sum(~isnan(x(:,ic)) & ~isnan(y(:,ic)));
        [spearRho(1,ic), pval(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Spearman', 'rows', 'complete');
        [rx(1,ic), pval_r(1,ic)] = corr(x(:,ic), y(:,ic), 'type', 'Pearson', 'rows', 'complete');
        [r(1,ic),b1(1,ic),b0(1,ic)] = regression(x(:,ic), y(:,ic), 'one');
        line([min(min(x))  max(max(x))], [b0(1,ic)+b1(1,ic)*min(min(x)) b0(1,ic)+b1(1,ic)*max(max(x))], 'color', col_lines(ic,:), 'linewidth', lw, 'linestyle', ls);
        xlim(xl); ylim(yl);
        xticks = -2:0.2:2; set(gca, 'XTick', xticks)
        yticks = -2:0.2:2;  set(gca, 'YTick', yticks);
    end
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


fh.Renderer = 'painters'; 
saveas(fh0, fullfile('dfi_experiment_figures', 'Paper_figures', 'iAF', ...
    'iAF_within', 'sd_params', 'Consistency_freq_sdparams_lcmv_svg.svg'))
close all



%% Bayes factor figures

% Fill in by hand!!!
BFs_dp = [0.1942, 0.1934, 0.4123];
BFs_c = [0.1792, 0.2646, 0.1771];

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
saveas(fh1, fullfile('dfi_experiment_figures', 'Paper_figures', 'iAF', ...
    'iAF_within', 'sd_params', 'Consistency_freq_sdparams_svg_lcmv_BFs.svg'))
close all


% eof

