% Compute dprime and bias for different beep contexts
% in yes-no and yes-no threshold task. 
% Creates Figure 1c (and corresponding statistics).
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
%       Stimulus overview
%
%     STIM1   Stim2      ID
%       V1    V1V2   =   2            V 
%      V1V2    V1    =   3            V 
%       A1    A1A2   =   4            A only
%      V1A1   V2A1   =   5            Fus
%      V2A1   V1A1   =   6            Fus 
%      A1A2    A1    =   7            A 
%      V1A2   V2A2   =   8            Fis
%      V2A2   V1A2   =   9            Fis 
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
% Last modified 02/05/2021


restoredefaultpath
clear all
close all

addpath(genpath('dfi'));
dfi_startup

data_dir = fullfile('dfi_experiment_data', 'data', 'experiment');
fig_dir = fullfile('dfi_experiment_figures', 'Paper_figures', 'iAF', 'beh');
mkdir(fig_dir)


%% -------2IFC--------

%% Statistics on d' and C data (2IFC task)

load(fullfile( data_dir, '2IFC', 'SD_params.mat' ));

dP_yesno = dP_adj;
C_yesno = C_adj;

dP_avg = squeeze(nanmean(dP_yesno(3:6,:,:)));
C_avg = squeeze(nanmean(C_yesno(3:6,:,:)));


% Create figure and save for inkscape: dPrime
yl = [0.9, 2.4];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = dP_avg';
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);

% get Cousineau within subject SE for plotting:

 
data = acc_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(acc_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(acc_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(acc_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(acc_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(acc_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(acc_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
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

%export_fig(fh, fullfile(fig_dir, fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all



% Create figure and save for inkscape: Criterion
yl = [-1.25, 1.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = C_avg';
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);

% get Cousineau within subject SE for plotting:
data = acc_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(acc_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(acc_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(acc_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(acc_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(acc_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(acc_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
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
%export_fig(fh, fullfile(fig_dir, fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all



% Create figure and save for inkscape: actual Criterion (before was bias_centre)
yl = [-1.25, 1.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = C_avg' + dP_avg'./2; %bias_centre + dp/2 = C
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);

% get Cousineau within subject SE for plotting:
data = acc_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(acc_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(acc_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(acc_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(acc_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(acc_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(acc_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
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
%export_fig(fh, fullfile(fig_dir, fname), '-tiff', '-m2.5');
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
mat_anova = dP_mat;
table_anova = array2table(mat_anova, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(table_anova,'Flash_fusion-Fission_illusion ~ 1','WithinDesign', ...
    withinDesign);

% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(mat_anova(:,2), mat_anova(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( mat_anova(:,1) - mat_anova(:,2) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,2)).^2)/2);
cohens_d(2) = nanmean( mat_anova(:,1) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d(3) = nanmean( mat_anova(:,2) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,2)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(mat_anova(:,2), mat_anova(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = mat_anova - repmat(nanmean(mat_anova,1),[size(mat_anova,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(mat_anova, 2);
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
for i = 1:size(mat_anova, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(mat_anova, 2)
    alpha = 0.001; % as it's conservative let this be quite small
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
mat_anova = C_mat;
table_anova = array2table(mat_anova, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});

% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});

% create model object
rm = fitrm(table_anova,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);

% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(mat_anova(:,2), mat_anova(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( mat_anova(:,1) - mat_anova(:,2) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,2)).^2)/2);
cohens_d(2) = nanmean( mat_anova(:,1) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d(3) = nanmean( mat_anova(:,2) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,2)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(mat_anova(:,2), mat_anova(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = mat_anova - repmat(nanmean(mat_anova,1),[size(mat_anova,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(mat_anova, 2);
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
for i = 1:size(mat_anova, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(mat_anova, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


%% stats Cnoise:

% Save data for SPSS
oldwd = cd;
C_mat = C_avg' + dP_avg'./2;
cd(data_dir, 'twoIFC_beh_stats'));
dlmwrite('Cnoise_mat_stats.txt', C_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
mat_anova = C_avg' + dP_avg'./2;
table_anova = array2table(mat_anova, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});

% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});

% create model object
rm = fitrm(table_anova,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);

% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(mat_anova(:,2), mat_anova(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( mat_anova(:,1) - mat_anova(:,2) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,2)).^2)/2);
cohens_d(2) = nanmean( mat_anova(:,1) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d(3) = nanmean( mat_anova(:,2) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,2)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(mat_anova(:,2), mat_anova(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = mat_anova - repmat(nanmean(mat_anova,1),[size(mat_anova,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(mat_anova, 2);
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
for i = 1:size(mat_anova, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(mat_anova, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


%% -------Yes-no intermediate SOAs--------

%% Statistics on d' and C data

load(fullfile( data_dir, 'yesno', 'SD_params.mat' ));

dP_yesno = dP_adj;
C_yesno = C_adj;

dP_avg = squeeze(nanmean(dP_yesno(3:6,:,:)));
C_avg = squeeze(nanmean(C_yesno(3:6,:,:)));


% Create figure and save for inkscape: dPrime
yl = [0.9, 2.4];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = dP_avg';
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);

% get Cousineau within subject SE for plotting:
data = acc_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(acc_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(acc_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(acc_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(acc_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(acc_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(acc_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
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
%export_fig(fh, fullfile(fig_dir, fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all


% Create figure and save for inkscape: Bias_centre
yl = [-1.25, 1.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = C_avg';
acc_mat(:,3) = -acc_mat(:,3);
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);

% get Cousineau within subject SE for plotting:
data = acc_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(acc_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(acc_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(acc_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(acc_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(acc_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(acc_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
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
%export_fig(fh, fullfile(fig_dir, fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all


% Create figure and save for inkscape: Criterion
yl = [-2.25, 2.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = C_avg' + dP_avg'./2; %bias_centre + dp/2 = C
acc_mat(:,3) = -C_avg(3,:)' + dP_avg(3,:)'./2;
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);

% get Cousineau within subject SE for plotting:
data = acc_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(acc_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(acc_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(acc_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(acc_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(acc_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(acc_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
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
%export_fig(fh, fullfile(fig_dir, fname), '-tiff', '-m2.5');
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
mat_anova = dP_mat;
table_anova = array2table(mat_anova, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});

% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});

% create model object
rm = fitrm(table_anova,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);

% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(mat_anova(:,2), mat_anova(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( mat_anova(:,1) - mat_anova(:,2) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,2)).^2)/2);
cohens_d(2) = nanmean( mat_anova(:,1) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d(3) = nanmean( mat_anova(:,2) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,2)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(mat_anova(:,2), mat_anova(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = mat_anova - repmat(nanmean(mat_anova,1),[size(mat_anova,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(mat_anova, 2);
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
for i = 1:size(mat_anova, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(mat_anova, 2);
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
C_mat(:,3) = -C_mat(:,3);
mat_anova = C_mat;
table_anova = array2table(mat_anova, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});

% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});

% create model object
rm = fitrm(table_anova,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);

% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(mat_anova(:,2), mat_anova(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( mat_anova(:,1) - mat_anova(:,2) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,2)).^2)/2);
cohens_d(2) = nanmean( mat_anova(:,1) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d(3) = nanmean( mat_anova(:,2) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,2)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(mat_anova(:,2), mat_anova(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = mat_anova - repmat(nanmean(mat_anova,1),[size(mat_anova,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(mat_anova, 2);
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
for i = 1:size(mat_anova, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(mat_anova, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


%% stats Cnoise:

% Save data for SPSS
oldwd = cd;
acc_mat = C_avg' + dP_avg'./2; %bias_centre + dp/2 = C
acc_mat(:,3) = -C_avg(3,:)' + dP_avg(3,:)'./2;
C_mat = acc_mat;
cd(fullfile(data_dir, 'yn_beh_stats'));
dlmwrite('Cnoise_mat_stats.txt', C_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
mat_anova = C_mat;
table_anova = array2table(mat_anova, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});

% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});

% create model object
rm = fitrm(table_anova,'Flash_fusion-Fission_illusion ~ 1','WithinDesign', withinDesign);

% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(mat_anova(:,2), mat_anova(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( mat_anova(:,1) - mat_anova(:,2) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,2)).^2)/2);
cohens_d(2) = nanmean( mat_anova(:,1) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d(3) = nanmean( mat_anova(:,2) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,2)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(mat_anova(:,2), mat_anova(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = mat_anova - repmat(nanmean(mat_anova,1),[size(mat_anova,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(mat_anova, 2);
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
for i = 1:size(mat_anova, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(mat_anova, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


%% -------Yes-no threshold-------

% Name of directory to save things to
an_fold = 'erps_ynt';

% experiment data folder
data_dir = fullfile('dfi_experiment_data', 'eeg_data', 'experiment');
fig_dir  = 'dfi_experiment_figures';

% useful variables
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
    '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
N        = 20;
task     = 'yn_threshold';

condvect = [3,6,8];
cond_labels = {'2F', 'Fus', 'Fis'};


% Load in trial indices to keep
load(fullfile(data_dir, 'eeg_stim_lock', 'ynt_trial_indeces_hands_matched_rt500.mat'))


% SOA distribution
tabulate(dall_stim.soa)
soavect = nan(20,3);
trlids = [3,6,8];
partvect = unique(dall_stim.partid);
for isubj = 1:20
    for icond = 1:3
        soavect(isubj,icond) = nanmean(unique(dall_stim.soa(dall_stim.partid == partvect(isubj) & dall_stim.trlid == trlids(icond))));
    end
end




%% Yes-no threshold: d', C

dP_mat = nan(20,3);
C_mat = nan(20,3);
pC_mat = nan(20,3);
actualC = nan(20,3);
for isubj = 1:20
    
    beh = dall_stim(dall_stim.partid == partvect(isubj),:);
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
acc_mat = dP_mat;
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);

% get Cousineau within subject SE for plotting:
data = acc_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(acc_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(acc_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(acc_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(acc_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(acc_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(acc_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
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
%export_fig(fh, fullfile(fig_dir, fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all


% Create figure and save for inkscape: Criterion (bias centre)
yl = [-1.25, 1.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = C_mat;
acc_mat(:,3) = -acc_mat(:,3);
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);

% get Cousineau within subject SE for plotting:
data = acc_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(acc_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(acc_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(acc_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(acc_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(acc_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(acc_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
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
%export_fig(fh, fullfile(fig_dir, fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all


% Create figure and save for inkscape: Criterion
yl = [-2.25, 2.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = actualC;
acc_mat(:,3) = -C_mat(:,3) + dP_mat(:,3)./2;
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);

% get Cousineau within subject SE for plotting:
data = acc_mat;
subj_mean = nanmean(data,2);
grand_mean = nanmean(subj_mean,1);
data_corr = data - repmat(subj_mean, [1,3]) + repmat(grand_mean, [size(data,1),size(data,2)]);
SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
SE_b = nanstd(data,0,1)./sqrt(sum(~isnan(data),1));

plot(1, nanmean(acc_mat(:,1)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,1)); hold on
plot(2, nanmean(acc_mat(:,2)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,2)); hold on
plot(3, nanmean(acc_mat(:,3)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', col_vect(:,3)); hold on
errorbar(1, nanmean(acc_mat(:,1)), SE_w(1), SE_w(1), 'Color', col_vect(:,1), 'LineWidth', 1.5);
errorbar(2, nanmean(acc_mat(:,2)), SE_w(2), SE_w(2), 'Color', col_vect(:,2), 'LineWidth', 1.5);
errorbar(3, nanmean(acc_mat(:,3)), SE_w(3), SE_w(3), 'Color', col_vect(:,3), 'LineWidth', 1.5);
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
%export_fig(fh, fullfile(fig_dir, fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile(fig_dir, [fname, '_svg.svg']))

close all





%% stats dP:

% RM-Anova
mat_anova = dP_mat;
table_anova = array2table(mat_anova, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});

% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});

% create model object
rm = fitrm(table_anova,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);

% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(mat_anova(:,2), mat_anova(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( mat_anova(:,1) - mat_anova(:,2) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,2)).^2)/2);
cohens_d(2) = nanmean( mat_anova(:,1) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d(3) = nanmean( mat_anova(:,2) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,2)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d
% sign rank test
[p,h,stats] = signrank(mat_anova(:,2), mat_anova(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = mat_anova - repmat(nanmean(mat_anova,1),[size(mat_anova,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(mat_anova, 2);
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
for i = 1:size(mat_anova, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(mat_anova, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


%% stats C:

% RM-Anova
C_mat(:,3) = -C_mat(:,3);
mat_anova = C_mat;
table_anova = array2table(mat_anova, 'VariableNames', ...
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});

% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});

% create model object
rm = fitrm(table_anova,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);

% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(mat_anova(:,2), mat_anova(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( mat_anova(:,1) - mat_anova(:,2) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,2)).^2)/2);
cohens_d(2) = nanmean( mat_anova(:,1) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d(3) = nanmean( mat_anova(:,2) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,2)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d
% sign rank test
[p,h,stats] = signrank(mat_anova(:,2), mat_anova(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = mat_anova - repmat(nanmean(mat_anova,1),[size(mat_anova,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(mat_anova, 2);
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
for i = 1:size(mat_anova, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(mat_anova, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


%% stats Cnoise:

% RM-Anova
mat_anova = acc_mat;
table_anova = array2table(mat_anova, 'VariableNames', ... 
    {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});

% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});

% create model object
rm = fitrm(table_anova,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);

% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(mat_anova(:,1), mat_anova(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(mat_anova(:,2), mat_anova(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( mat_anova(:,1) - mat_anova(:,2) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,2)).^2)/2);
cohens_d(2) = nanmean( mat_anova(:,1) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,1)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d(3) = nanmean( mat_anova(:,2) - mat_anova(:,3) ) / ...
    sqrt((nanstd(mat_anova(:,2)).^2 + nanstd(mat_anova(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(mat_anova(:,2), mat_anova(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = mat_anova - repmat(nanmean(mat_anova,1),[size(mat_anova,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(mat_anova, 2);
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
for i = 1:size(mat_anova, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(mat_anova, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


% eof

