% Behavioural statistics and figures for pre-stim manuscript (SIFI project)
%
% Mostly recycled from post-stim and psychophysics projects
%
% sb, sbuergers@gmail.com, 05/12/2019


clear all
close all


% access matlab file exchange functions
addpath(genpath('D:\matlab_file_exchange'));

data_dir = 'D:\dfi_experiment_data\data\experiment';



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
%distributionPlot_DM([acc_mat],'globalNorm',true,'colormap',1-gray(128),'histOpt',1, 'showMM', 0,'addSpread',2, 'distWidth', 0.95); hold on % histOpt=2 works better for uniform distributions than the default
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
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
mkdir('D:\dfi_experiment_figures\Paper_figures\iAF\beh');
export_fig(fh, fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', [fname, '_svg.svg']))

close all



% Create figure and save for inkscape: Criterion
yl = [-1.25, 1.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = C_avg';
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
%distributionPlot_DM([acc_mat],'globalNorm',true,'colormap',1-gray(128),'histOpt',1, 'showMM', 0,'addSpread',2, 'distWidth', 0.95); hold on % histOpt=2 works better for uniform distributions than the default
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
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
mkdir('D:\dfi_experiment_figures\Paper_figures\iAF\beh');
export_fig(fh, fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', [fname, '_svg.svg']))

close all



% Create figure and save for inkscape: actual Criterion (before was bias_centre)
yl = [-1.25, 1.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = C_avg' + dP_avg'./2; %bias_centre + dp/2 = C
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
%distributionPlot_DM([acc_mat],'globalNorm',true,'colormap',1-gray(128),'histOpt',1, 'showMM', 0,'addSpread',2, 'distWidth', 0.95); hold on % histOpt=2 works better for uniform distributions than the default
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
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
mkdir('D:\dfi_experiment_figures\Paper_figures\iAF\beh');
export_fig(fh, fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', [fname, '_svg.svg']))

close all





%% stats dP:

% Save data for SPSS
oldwd = cd;
dP_mat = dP_avg';
mkdir('D:\dfi_experiment_data\data\experiment\twoIFC_beh_stats')
cd('D:\dfi_experiment_data\data\experiment\twoIFC_beh_stats');
dlmwrite('dP_mat_stats.txt', dP_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
threshold_matrix = dP_mat;
absbiasTable = array2table(threshold_matrix, 'VariableNames', {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(absbiasTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');
% Note that given that this table is balanced there is no difference
% between type I, II or III sum of squares, this only occurs when data is
% unbalanced. Generally II and III are preferable over I, as they take into
% account other main effects assuming no interactions and other main
% effects and interactions respectively (if there is an interaction main
% effects are not really interpretable). See:
% https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(threshold_matrix(:,2), threshold_matrix(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,2) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,2)).^2)/2);
cohens_d(2) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d(3) = nanmean( threshold_matrix(:,2) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,2)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(threshold_matrix(:,2), threshold_matrix(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = threshold_matrix - repmat(nanmean(threshold_matrix,1),[size(threshold_matrix,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(threshold_matrix, 2);
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
for i = 1:size(threshold_matrix, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(threshold_matrix, 2);
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
cd('D:\dfi_experiment_data\data\experiment\twoIFC_beh_stats');
dlmwrite('C_mat_stats.txt', C_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
threshold_matrix = C_mat;
absbiasTable = array2table(threshold_matrix, 'VariableNames', {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(absbiasTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');
% Note that given that this table is balanced there is no difference
% between type I, II or III sum of squares, this only occurs when data is
% unbalanced. Generally II and III are preferable over I, as they take into
% account other main effects assuming no interactions and other main
% effects and interactions respectively (if there is an interaction main
% effects are not really interpretable). See:
% https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(threshold_matrix(:,2), threshold_matrix(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,2) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,2)).^2)/2);
cohens_d(2) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d(3) = nanmean( threshold_matrix(:,2) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,2)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(threshold_matrix(:,2), threshold_matrix(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = threshold_matrix - repmat(nanmean(threshold_matrix,1),[size(threshold_matrix,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(threshold_matrix, 2);
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
for i = 1:size(threshold_matrix, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(threshold_matrix, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


%% stats Cnoise:

% Save data for SPSS
oldwd = cd;
C_mat = C_avg' + dP_avg'./2;
cd('D:\dfi_experiment_data\data\experiment\twoIFC_beh_stats');
dlmwrite('Cnoise_mat_stats.txt', C_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
threshold_matrix = C_avg' + dP_avg'./2;
absbiasTable = array2table(threshold_matrix, 'VariableNames', {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(absbiasTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');
% Note that given that this table is balanced there is no difference
% between type I, II or III sum of squares, this only occurs when data is
% unbalanced. Generally II and III are preferable over I, as they take into
% account other main effects assuming no interactions and other main
% effects and interactions respectively (if there is an interaction main
% effects are not really interpretable). See:
% https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(threshold_matrix(:,2), threshold_matrix(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,2) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,2)).^2)/2);
cohens_d(2) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d(3) = nanmean( threshold_matrix(:,2) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,2)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(threshold_matrix(:,2), threshold_matrix(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = threshold_matrix - repmat(nanmean(threshold_matrix,1),[size(threshold_matrix,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(threshold_matrix, 2);
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
for i = 1:size(threshold_matrix, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(threshold_matrix, 2);
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
acc_mat = dP_avg';
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
%distributionPlot_DM([acc_mat],'globalNorm',true,'colormap',1-gray(128),'histOpt',1, 'showMM', 0,'addSpread',2, 'distWidth', 0.95); hold on % histOpt=2 works better for uniform distributions than the default
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
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
mkdir('D:\dfi_experiment_figures\Paper_figures\iAF\beh');
export_fig(fh, fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', [fname, '_svg.svg']))

close all



% Create figure and save for inkscape: Criterion
yl = [-1.25, 1.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = C_avg';
acc_mat(:,3) = -acc_mat(:,3);
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
%distributionPlot_DM([acc_mat],'globalNorm',true,'colormap',1-gray(128),'histOpt',1, 'showMM', 0,'addSpread',2, 'distWidth', 0.95); hold on % histOpt=2 works better for uniform distributions than the default
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
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
mkdir('D:\dfi_experiment_figures\Paper_figures\iAF\beh');
export_fig(fh, fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', [fname, '_svg.svg']))

close all


% Create figure and save for inkscape: Criterion
yl = [-2.25, 2.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = C_avg' + dP_avg'./2; %bias_centre + dp/2 = C
acc_mat(:,3) = -C_avg(3,:)' + dP_avg(3,:)'./2;
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
%distributionPlot_DM([acc_mat],'globalNorm',true,'colormap',1-gray(128),'histOpt',1, 'showMM', 0,'addSpread',2, 'distWidth', 0.95); hold on % histOpt=2 works better for uniform distributions than the default
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
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
mkdir('D:\dfi_experiment_figures\Paper_figures\iAF\beh');
export_fig(fh, fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', [fname, '_svg.svg']))

close all

%% stats dP:

% Save data for SPSS
oldwd = cd;
dP_mat = dP_avg';
mkdir('D:\dfi_experiment_data\data\experiment\yn_beh_stats')
cd('D:\dfi_experiment_data\data\experiment\yn_beh_stats');
dlmwrite('dP_mat_stats.txt', dP_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
threshold_matrix = dP_mat;
absbiasTable = array2table(threshold_matrix, 'VariableNames', {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(absbiasTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');
% Note that given that this table is balanced there is no difference
% between type I, II or III sum of squares, this only occurs when data is
% unbalanced. Generally II and III are preferable over I, as they take into
% account other main effects assuming no interactions and other main
% effects and interactions respectively (if there is an interaction main
% effects are not really interpretable). See:
% https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(threshold_matrix(:,2), threshold_matrix(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,2) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,2)).^2)/2);
cohens_d(2) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d(3) = nanmean( threshold_matrix(:,2) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,2)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(threshold_matrix(:,2), threshold_matrix(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = threshold_matrix - repmat(nanmean(threshold_matrix,1),[size(threshold_matrix,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(threshold_matrix, 2);
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
for i = 1:size(threshold_matrix, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(threshold_matrix, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output

%% stats C:

% Save data for SPSS
oldwd = cd;
C_mat = C_avg';
cd('D:\dfi_experiment_data\data\experiment\yn_beh_stats');
dlmwrite('C_mat_stats.txt', C_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
C_mat(:,3) = -C_mat(:,3);
threshold_matrix = C_mat;
absbiasTable = array2table(threshold_matrix, 'VariableNames', {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(absbiasTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');
% Note that given that this table is balanced there is no difference
% between type I, II or III sum of squares, this only occurs when data is
% unbalanced. Generally II and III are preferable over I, as they take into
% account other main effects assuming no interactions and other main
% effects and interactions respectively (if there is an interaction main
% effects are not really interpretable). See:
% https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(threshold_matrix(:,2), threshold_matrix(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,2) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,2)).^2)/2);
cohens_d(2) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d(3) = nanmean( threshold_matrix(:,2) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,2)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(threshold_matrix(:,2), threshold_matrix(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = threshold_matrix - repmat(nanmean(threshold_matrix,1),[size(threshold_matrix,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(threshold_matrix, 2);
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
for i = 1:size(threshold_matrix, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(threshold_matrix, 2);
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
cd('D:\dfi_experiment_data\data\experiment\yn_beh_stats');
dlmwrite('Cnoise_mat_stats.txt', C_mat, 'delimiter', '\t');
cd(oldwd)

% RM-Anova
threshold_matrix = C_mat;
absbiasTable = array2table(threshold_matrix, 'VariableNames', {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(absbiasTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');
% Note that given that this table is balanced there is no difference
% between type I, II or III sum of squares, this only occurs when data is
% unbalanced. Generally II and III are preferable over I, as they take into
% account other main effects assuming no interactions and other main
% effects and interactions respectively (if there is an interaction main
% effects are not really interpretable). See:
% https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(threshold_matrix(:,2), threshold_matrix(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,2) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,2)).^2)/2);
cohens_d(2) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d(3) = nanmean( threshold_matrix(:,2) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,2)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d

% sign rank test
[p,h,stats] = signrank(threshold_matrix(:,2), threshold_matrix(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = threshold_matrix - repmat(nanmean(threshold_matrix,1),[size(threshold_matrix,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(threshold_matrix, 2);
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
for i = 1:size(threshold_matrix, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(threshold_matrix, 2);
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
load(fullfile('dfi_experiment_data', 'data', 'experiment', 'eeg_stim_lock', 'ynt_trial_indeces_hands_matched_rt500.mat'))


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
%distributionPlot_DM([acc_mat],'globalNorm',true,'colormap',1-gray(128),'histOpt',1, 'showMM', 0,'addSpread',2, 'distWidth', 0.95); hold on % histOpt=2 works better for uniform distributions than the default
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
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
mkdir('D:\dfi_experiment_figures\Paper_figures\iAF\beh');
export_fig(fh, fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', [fname, '_svg.svg']))

close all



% Create figure and save for inkscape: Criterion
yl = [-1.25, 1.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = C_mat;
acc_mat(:,3) = -acc_mat(:,3);
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
%distributionPlot_DM([acc_mat],'globalNorm',true,'colormap',1-gray(128),'histOpt',1, 'showMM', 0,'addSpread',2, 'distWidth', 0.95); hold on % histOpt=2 works better for uniform distributions than the default
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
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
mkdir('D:\dfi_experiment_figures\Paper_figures\iAF\beh');
export_fig(fh, fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', [fname, '_svg.svg']))

close all




% Create figure and save for inkscape: Criterion
yl = [-2.25, 2.25];
cond_labels = {'2F', 'Fus', 'Fis'};
col_vect = [0 0.6 0; 0 0 1; 1 0 0]';
acc_mat = actualC;
acc_mat(:,3) = -C_mat(:,3) + dP_mat(:,3)./2;
fh = figure('color', [1 1 1], 'Position', [0, 0, 140, 140]);
%distributionPlot_DM([acc_mat],'globalNorm',true,'colormap',1-gray(128),'histOpt',1, 'showMM', 0,'addSpread',2, 'distWidth', 0.95); hold on % histOpt=2 works better for uniform distributions than the default
% get Cousineau within subject SE for plotting:
% Cancel out between subject variability by subtracting the subject
% mean from each subject and then adding the grand mean
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
mkdir('D:\dfi_experiment_figures\Paper_figures\iAF\beh');
export_fig(fh, fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', fname), '-tiff', '-m2.5');
fh.Renderer = 'painters'; saveas(fh,  fullfile('D:\dfi_experiment_figures\Paper_figures\iAF\beh', [fname, '_svg.svg']))

close all





%% stats dP:

% RM-Anova
threshold_matrix = dP_mat;
absbiasTable = array2table(threshold_matrix, 'VariableNames', {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(absbiasTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');
% Note that given that this table is balanced there is no difference
% between type I, II or III sum of squares, this only occurs when data is
% unbalanced. Generally II and III are preferable over I, as they take into
% account other main effects assuming no interactions and other main
% effects and interactions respectively (if there is an interaction main
% effects are not really interpretable). See:
% https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(threshold_matrix(:,2), threshold_matrix(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,2) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,2)).^2)/2);
cohens_d(2) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d(3) = nanmean( threshold_matrix(:,2) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,2)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d
% sign rank test
[p,h,stats] = signrank(threshold_matrix(:,2), threshold_matrix(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = threshold_matrix - repmat(nanmean(threshold_matrix,1),[size(threshold_matrix,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(threshold_matrix, 2);
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
for i = 1:size(threshold_matrix, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(threshold_matrix, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output




%% stats C:

% RM-Anova
C_mat(:,3) = -C_mat(:,3);
threshold_matrix = C_mat;
absbiasTable = array2table(threshold_matrix, 'VariableNames', {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(absbiasTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');
% Note that given that this table is balanced there is no difference
% between type I, II or III sum of squares, this only occurs when data is
% unbalanced. Generally II and III are preferable over I, as they take into
% account other main effects assuming no interactions and other main
% effects and interactions respectively (if there is an interaction main
% effects are not really interpretable). See:
% https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(threshold_matrix(:,2), threshold_matrix(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,2) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,2)).^2)/2);
cohens_d(2) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d(3) = nanmean( threshold_matrix(:,2) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,2)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d
% sign rank test
[p,h,stats] = signrank(threshold_matrix(:,2), threshold_matrix(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = threshold_matrix - repmat(nanmean(threshold_matrix,1),[size(threshold_matrix,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(threshold_matrix, 2);
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
for i = 1:size(threshold_matrix, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(threshold_matrix, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output



%% stats Cnoise:

% RM-Anova
threshold_matrix = acc_mat;
absbiasTable = array2table(threshold_matrix, 'VariableNames', {'Flash_fusion', 'Fusion_illusion', 'Fission_illusion'});
% Create design table
withinDesign = table(...
    categorical({'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'}),...
    'VariableNames',{'illusion'});
% create model object
rm = fitrm(absbiasTable,'Flash_fusion-Fission_illusion ~ 1','WithinDesign',withinDesign);
% run repeated measures anova (this gives the same results as ezANOVA with
% type 2 sum of squares in R!
ranovaTable = ranova(rm,'WithinModel','illusion');
% Note that given that this table is balanced there is no difference
% between type I, II or III sum of squares, this only occurs when data is
% unbalanced. Generally II and III are preferable over I, as they take into
% account other main effects assuming no interactions and other main
% effects and interactions respectively (if there is an interaction main
% effects are not really interpretable). See:
% https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/

% Test sphericity assumption.
fprintf('\nMauchly`s test of sphericity:\n')
tbl = mauchly(rm)

% t-tests
fprintf('\ncompare flash fusion to fusion illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,2))
fprintf('\ncompare flash fusion to fission illusion\n')
[H,P,CI,STATS] = ttest(threshold_matrix(:,1), threshold_matrix(:,3))
fprintf('\ncompare fusion illusion to fission illusion\n')
[~,P,CI,STATS] = ttest(threshold_matrix(:,2), threshold_matrix(:,3))

% effect size
cohens_d = nan(3,1);
cohens_d(1) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,2) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,2)).^2)/2);
cohens_d(2) = nanmean( threshold_matrix(:,1) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,1)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d(3) = nanmean( threshold_matrix(:,2) - threshold_matrix(:,3) ) / ...
    sqrt((nanstd(threshold_matrix(:,2)).^2 + nanstd(threshold_matrix(:,3)).^2)/2);
cohens_d
% sign rank test
[p,h,stats] = signrank(threshold_matrix(:,2), threshold_matrix(:,3));

% test assumptions, plot residuals!
% 1.) Homoscedasticity
figure('color', 'w')
residuals = threshold_matrix - repmat(nanmean(threshold_matrix,1),[size(threshold_matrix,1),1]);
titles = {'Flash_fusion'; 'Fusion_illusion'; 'Fission_illusion'};
for i = 1:size(threshold_matrix, 2);
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
for i = 1:size(threshold_matrix, 2);
    subplot(1,3,i);
    qqplot(residuals(:,i));
    title(titles{i});
    if i <= 4, xlabel(''); end;
    if ~ismember(i, [1,5,9,13]), ylabel(''); end;
end
suptitle('Residual plots');

% Test normality with shapiro wilk test (very conservative test!)
[H, pValue, SWstatistic] = deal(nan(size(residuals,2),1));
for i = 1:size(threshold_matrix, 2);
    alpha = 0.001; %as it's conservative let this be quite small
    [H(i), pValue(i), SWstatistic(i)] = swtest(residuals(:,i), alpha);
end
shapiro_output = [H, pValue, SWstatistic];
shapiro_output


