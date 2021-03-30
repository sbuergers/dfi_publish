% Behavioural statistics and figures for dprime and bias.
% "Creates" Supplementary Table 3
%
% Parent script(s): 
%   dfi_sdt_analyses.m
%   
% Children script(s): 
%   None
%
% Sibling script(s):
%   None
%
%
% DETAILS
%
% Note that comparing different auditory contexts is problematic for
% the yes-no threshold task, because SOAs were titrated separately
% for them and may not be matched. 
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
% Last modified 30/03/2021


restoredefaultpath
clear all
close all


addpath(genpath('dfi'))
dfi_startup
data_dir = fullfile('dfi_experiment_data', 'data', 'experiment');
fig_dir = fullfile('dfi_experiment_figures', 'Paper_figures', 'iAF', 'beh');
mkdir(fig_dir);


%% -------2IFC--------

%% Statistics on d' and C data (2IFC task)

load(fullfile( data_dir, '2IFC', 'SD_params.mat' ));

dP_2ifc = dP_adj;
C_2ifc = C_adj;

dP_avg = squeeze(nanmean(dP_2ifc(3:6,:,:)));
C_avg = squeeze(nanmean(C_2ifc(3:6,:,:)));


% Create figure and save for inkscape: dPrime
yl = [0.9, 2.4];
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
x_mat = dP_avg';
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
data = x_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(x_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(x_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(x_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(x_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(x_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(x_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
set(gca, 'xticklabel', [])
%set(gca, 'yticklabel', [])
box off
set(gca,'TickDir','out')
yticks = 0:0.5:100;  set(gca, 'YTick', yticks);
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')
ylim(yl)
xlim([0.5 3.5])

% Export figures
fname = 'beh_dP_adj_2ifc';
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all



% Create figure and save for inkscape: Criterion
yl = [-1.25, 1.25];
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
x_mat = C_avg';
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
data = x_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(x_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(x_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(x_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(x_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(x_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(x_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
set(gca, 'xticklabel', [])
%set(gca, 'yticklabel', [])
box off
set(gca,'TickDir','out')
yticks = -100:0.5:100;  set(gca, 'YTick', yticks);
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')
ylim(yl)
xlim([0.5 3.5])

% Export figures
fname = 'beh_C_adj_2ifc';
mkdir(fig_dir);
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all



% Create figure and save for inkscape: actual Criterion (before was bias_centre)
yl = [-1.25, 1.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
x_mat = C_avg' + dP_avg'./2; %bias_centre + dp/2 = C
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
data = x_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(x_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(x_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(x_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(x_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(x_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(x_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
set(gca, 'xticklabel', [])
%set(gca, 'yticklabel', [])
box off
set(gca,'TickDir','out')
yticks = -100:0.5:100;  set(gca, 'YTick', yticks);
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')
ylim(yl)
xlim([0.5 3.5])

% Export figures
fname = 'beh_actualC_adj_2ifc';
mkdir(fig_dir);
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all





%% stats dP:

% Save data for SPSS
oldwd = cd;
dP_mat = dP_avg';
mkdir(fullfile(data_dir, 'twoIFC_beh_stats'))
cd(fullfile(data_dir, 'twoIFC_beh_stats'));
dlmwrite('dP_mat_stats.txt', dP_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
x = dP_mat;
xTable = array2table(x, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(xTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(x(:,2), x(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( x(:,1) - x(:,2) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,2)).^2)/2);
cohens_d(2) = nanmean( x(:,1) - x(:,3) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d(3) = nanmean( x(:,2) - x(:,3) ) / ...
    sqrt((nanstd(x(:,2)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(x(:,2), x(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = x - repmat(nanmean(x,1),[size(x,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(x, 2);
    subplot(1,3,i);
    scatter(1:size(residuals,1), residuals(:,i), 'ro', 'filled');
    title(titles{i});
    if i > 4, xlabel('Subject'); end;
    if ismember(i, [1,5,9,13]), ylabel('Residual value'); end;
end
suptitle('Residual plots');

% Plot residual quantiles versus normal distribution
% 2.) Normality
figure('color', 'w');
for i = 1:size(x, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(x, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


%% stats C (disregard SOA, disregard Context):
% i.e. simply check if the average Interval bias for a given sensory
% context is significantly different from zero
for icond = 1:3
    d = C_avg(icond,:);
    
    fprintf('\nCondition %i', icond);
    
    % ttest
    [H,P,CI,STATS] = ttest(d)
    
    % sign rank test
    [p,h,stats] = signrank(d)
end


%% stats C:

% Save data for SPSS
oldwd = cd;
C_mat = C_avg';
cd(fullfile(data_dir, 'twoIFC_beh_stats'));
dlmwrite('C_mat_stats.txt', C_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
x = C_mat;
xTable = array2table(x, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(xTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(x(:,2), x(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( x(:,1) - x(:,2) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,2)).^2)/2);
cohens_d(2) = nanmean( x(:,1) - x(:,3) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d(3) = nanmean( x(:,2) - x(:,3) ) / ...
    sqrt((nanstd(x(:,2)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(x(:,2), x(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = x - repmat(nanmean(x,1),[size(x,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(x, 2);
    subplot(1,3,i);
    scatter(1:size(residuals,1), residuals(:,i), 'ro', 'filled');
    title(titles{i});
    if i > 4, xlabel('Subject'); end;
    if ismember(i, [1,5,9,13]), ylabel('Residual value'); end;
end
suptitle('Residual plots');

% Plot residual quantiles versus normal distribution
% 2.) Normality
figure('color', 'w');
for i = 1:size(x, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(x, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


%% stats Cnoise:

% Save data for SPSS
oldwd = cd;
C_mat = C_avg' + dP_avg'./2;
cd(fullfile(data_dir, 'twoIFC_beh_stats'));
dlmwrite('Cnoise_mat_stats.txt', C_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
x = C_avg' + dP_avg'./2;
xTable = array2table(x, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(xTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(x(:,2), x(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( x(:,1) - x(:,2) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,2)).^2)/2);
cohens_d(2) = nanmean( x(:,1) - x(:,3) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d(3) = nanmean( x(:,2) - x(:,3) ) / ...
    sqrt((nanstd(x(:,2)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(x(:,2), x(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = x - repmat(nanmean(x,1),[size(x,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(x, 2);
    subplot(1,3,i);
    scatter(1:size(residuals,1), residuals(:,i), 'ro', 'filled');
    title(titles{i});
    if i > 4, xlabel('Subject'); end;
    if ismember(i, [1,5,9,13]), ylabel('Residual value'); end;
end
suptitle('Residual plots');

% Plot residual quantiles versus normal distribution
% 2.) Normality
figure('color', 'w');
for i = 1:size(x, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(x, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output



%% -------Yes-no intermediate SOAs--------

%% Statistics on d' and C data (not SVD outcomes)

load(fullfile( data_dir, 'yesno', 'SD_params.mat' ));

dP_yesno = dP_adj;
C_yesno = C_adj;

dP_avg = squeeze(nanmean(dP_yesno(3:6,:,:)));
C_avg = squeeze(nanmean(C_yesno(3:6,:,:)));


% Create figure and save for inkscape: dPrime
yl = [0.9, 2.4];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
x_mat = dP_avg';
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
data = x_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(x_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(x_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(x_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(x_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(x_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(x_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
set(gca, 'xticklabel', [])
%set(gca, 'yticklabel', [])
box off
set(gca,'TickDir','out')
yticks = 0:0.5:100;  set(gca, 'YTick', yticks);
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')
ylim(yl)
xlim([0.5 3.5])

% Export figures
fname = 'beh_dP_adj';
mkdir(fig_dir);
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all



% Create figure and save for inkscape: Criterion
yl = [-1.25, 1.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
x_mat = C_avg';
x_mat(:,3) = -x_mat(:,3);
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
%distributionPlot_DM([acc_mat],'globalNorm',true,'colormap',1-gray(128),'histOpt',1, 'showMM', 0,'addSpread',2, 'distWidth', 0.95); hold on % histOpt=2 works better for uniform distributions than the default
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
data = x_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(x_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(x_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(x_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(x_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(x_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(x_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
set(gca, 'xticklabel', [])
%set(gca, 'yticklabel', [])
box off
set(gca,'TickDir','out')
yticks = -100:0.5:100;  set(gca, 'YTick', yticks);
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')
ylim(yl)
xlim([0.5 3.5])

% Export figures
fname = 'beh_C_adj_2f2s_as_noise';
mkdir(fig_dir);
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all


% Create figure and save for inkscape: Criterion
yl = [-2.25, 2.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
x_mat = C_avg' + dP_avg'./2; %bias_centre + dp/2 = C
x_mat(:,3) = -C_avg(3,:)' + dP_avg(3,:)'./2;
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
%distributionPlot_DM([acc_mat],'globalNorm',true,'colormap',1-gray(128),'histOpt',1, 'showMM', 0,'addSpread',2, 'distWidth', 0.95); hold on % histOpt=2 works better for uniform distributions than the default
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
data = x_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(x_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(x_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(x_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(x_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(x_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(x_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
set(gca, 'xticklabel', [])
%set(gca, 'yticklabel', [])
box off
set(gca,'TickDir','out')
yticks = -100:0.5:100;  set(gca, 'YTick', yticks);
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')
ylim(yl)
xlim([0.5 3.5])

% Export figures
fname = 'beh_actualC_adj_2f2s_as_noise';
mkdir(fig_dir);
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all

%% stats dP:

% Save data for SPSS
oldwd = cd;
dP_mat = dP_avg';
mkdir(fullfile(data_dir, 'yn_beh_stats'))
cd(fullfile(data_dir, 'yn_beh_stats'));
dlmwrite('dP_mat_stats.txt', dP_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
x = dP_mat;
xTable = array2table(x, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(xTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(x(:,2), x(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( x(:,1) - x(:,2) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,2)).^2)/2);
cohens_d(2) = nanmean( x(:,1) - x(:,3) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d(3) = nanmean( x(:,2) - x(:,3) ) / ...
    sqrt((nanstd(x(:,2)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = x - repmat(nanmean(x,1),[size(x,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(x, 2);
    subplot(1,3,i);
    scatter(1:size(residuals,1), residuals(:,i), 'ro', 'filled');
    title(titles{i});
    if i > 4, xlabel('Subject'); end;
    if ismember(i, [1,5,9,13]), ylabel('Residual value'); end;
end
suptitle('Residual plots');

% Plot residual quantiles versus normal distribution
% 2.) Normality
figure('color', 'w');
for i = 1:size(x, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(x, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output

%% stats C:

% Save data for SPSS
oldwd = cd;
C_mat = C_avg';
cd(fullfile(data_dir, 'yn_beh_stats'));
dlmwrite('C_mat_stats.txt', C_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
%C_mat(:,3) = -C_mat(:,3);
x = C_mat;
xTable = array2table(x, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(xTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(x(:,2), x(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( x(:,1) - x(:,2) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,2)).^2)/2);
cohens_d(2) = nanmean( x(:,1) - x(:,3) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d(3) = nanmean( x(:,2) - x(:,3) ) / ...
    sqrt((nanstd(x(:,2)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(x(:,2), x(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = x - repmat(nanmean(x,1),[size(x,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(x, 2);
    subplot(1,3,i);
    scatter(1:size(residuals,1), residuals(:,i), 'ro', 'filled');
    title(titles{i});
    if i > 4, xlabel('Subject'); end;
    if ismember(i, [1,5,9,13]), ylabel('Residual value'); end;
end
suptitle('Residual plots');

% Plot residual quantiles versus normal distribution
% 2.) Normality
figure('color', 'w');
for i = 1:size(x, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(x, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output



%% -------Yes-no threshold-------

% Load in trial indices to keep
load(fullfile(data_dir, 'd701to727_ynt.mat'))


% delete no response and bad response trials
draw  = dall;
dall(dall.resp  == 0, :) = [];
dall(dall.badRT ~= 0, :) = [];
dall(dall.RT < 0.1,:)    = [];


% SOA distribution
tabulate(dall.soa)
soavect = nan(20,3);
trlids = [3,6,8];
partvect = unique(dall.partid);
for isubj = 1:20
    for icond = 1:3
        soavect(isubj,icond) = nanmean(unique(dall.soa(...
            dall.partid == partvect(isubj) & dall.trlid == trlids(icond))));
    end
end




%% Yes-no threshold: d', C

dP_mat = nan(20,3);
C_mat = nan(20,3);
pC_mat = nan(20,3);
actualC = nan(20,3);
for isubj = 1:20
    
    beh = dall(dall.partid == partvect(isubj),:);
    % 3.1) Type 1 signal detection parameters
    %S %N
    [~,sdm23] = dfi_SDM(beh,  3, 2, 0);
    [~,sdm56] = dfi_SDM(beh,  6, 5, 0);
    [~,sdm89] = dfi_SDM(beh,  9, 8, 0);
    
    dP_mat(isubj,1) = sdm23.dP_adj;
    C_mat(isubj,1) = sdm23.C_adj;
    pC_mat(isubj,1) = sdm23.pC_adj;
    
    dP_mat(isubj,2) = sdm56.dP_adj;
    C_mat(isubj,2) = sdm56.C_adj;
    pC_mat(isubj,2) = sdm56.pC_adj;
    
    dP_mat(isubj,3) = sdm89.dP_adj;
    C_mat(isubj,3) = sdm89.C_adj;
    pC_mat(isubj,3) = sdm89.pC_adj;
    
    actualC(isubj,1) = sdm23.actualC;
    actualC(isubj,2) = sdm56.actualC;
    actualC(isubj,3) = sdm89.actualC;
    
end

% Create figure and save for inkscape: dPrime
yl = [0.9, 2.4];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
x_mat = dP_mat;
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
data = x_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(x_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(x_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(x_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(x_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(x_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(x_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
set(gca, 'xticklabel', [])
%set(gca, 'yticklabel', [])
box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')
xlim([0.5 3.5])
ylim(yl)



% Export figures
fname = 'dP_adj_ynt';
mkdir(fig_dir);
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all



% Create figure and save for inkscape: Criterion
yl = [-1.25, 1.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
x_mat = C_mat;
x_mat(:,3) = -x_mat(:,3);
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
data = x_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(x_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(x_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(x_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(x_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(x_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(x_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
set(gca, 'xticklabel', [])
%set(gca, 'yticklabel', [])
box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')
xlim([0.5 3.5])
ylim(yl)

% Export figures
fname = 'C_adj_ynt_2f2s_as_noise';
mkdir(fig_dir);
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all




% Create figure and save for inkscape: Criterion
yl = [-2.25, 2.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
x_mat = actualC;
x_mat(:,3) = -C_mat(:,3) + dP_mat(:,3)./2;
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
data = x_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(x_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(x_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(x_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(x_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(x_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(x_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
set(gca, 'xticklabel', [])
%set(gca, 'yticklabel', [])
box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')
xlim([0.5 3.5])
ylim(yl)

% Export figures
fname = 'actualC_adj_ynt_2f2s_as_noise';
mkdir(fig_dir);
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all





%% stats dP:

% RM-Anova
x = dP_mat;
xTable = array2table(x, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(xTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(x(:,2), x(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( x(:,1) - x(:,2) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,2)).^2)/2);
cohens_d(2) = nanmean( x(:,1) - x(:,3) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d(3) = nanmean( x(:,2) - x(:,3) ) / ...
    sqrt((nanstd(x(:,2)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d
% sign rank test
[p,h,stats] = signrank(x(:,2), x(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = x - repmat(nanmean(x,1),[size(x,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(x, 2);
    subplot(1,3,i);
    scatter(1:size(residuals,1), residuals(:,i), 'ro', 'filled');
    title(titles{i});
    if i > 4, xlabel('Subject'); end;
    if ismember(i, [1,5,9,13]), ylabel('Residual value'); end;
end
suptitle('Residual plots');

% Plot residual quantiles versus normal distribution
% 2.) Normality
figure('color', 'w');
for i = 1:size(x, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(x, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output




%% stats C:

% RM-Anova
%C_mat(:,3) = -C_mat(:,3);
x = C_mat;
xTable = array2table(x, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(xTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(x(:,1), x(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(x(:,2), x(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( x(:,1) - x(:,2) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,2)).^2)/2);
cohens_d(2) = nanmean( x(:,1) - x(:,3) ) / ...
    sqrt((nanstd(x(:,1)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d(3) = nanmean( x(:,2) - x(:,3) ) / ...
    sqrt((nanstd(x(:,2)).^2 + nanstd(x(:,3)).^2)/2);
cohens_d
% sign rank test
[p,h,stats] = signrank(x(:,2), x(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = x - repmat(nanmean(x,1),[size(x,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(x, 2);
    subplot(1,3,i);
    scatter(1:size(residuals,1), residuals(:,i), 'ro', 'filled');
    title(titles{i});
    if i > 4, xlabel('Subject'); end;
    if ismember(i, [1,5,9,13]), ylabel('Residual value'); end;
end
suptitle('Residual plots');

% Plot residual quantiles versus normal distribution
% 2.) Normality
figure('color', 'w');
for i = 1:size(x, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(x, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


% eof

