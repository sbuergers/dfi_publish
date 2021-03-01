% Create figure of psychometric functions
%
% Parent script(s): 
%   dfi_beta_binom_combine_joined_and_sep_fits.m
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
% None
% 
% ---
% Steffen Buergers, sbuergers@gmail.com,
% Last modified Feb. 2021


clear all
close all


% access matlab file exchange functions
addpath(genpath('D:\matlab_file_exchange'));



%% *** Get data YN ***

% Paths, and data
cd('D:\')
addpath(genpath('dfi'));
data_dir = 'D:\dfi_experiment_data\data\experiment';
fig_dir  = 'D:\dfi_experiment_figures\group_figures\beh\yn\701to727';
load(fullfile(data_dir, 'd701to727_yn.mat'))
figdir = 'D:\dfi_experiment_figures\PFs\beta_binom_weibull';
dyn = dall; clear dall

savefigures = 1;


S2_conditions = [3,6,8,9];


%% *** Preprocessing and participant selection ***
% delete no response and bad response trials
draw_yn  = dyn;
dyn(dyn.resp  == 0, :) = [];
dyn(dyn.badRT ~= 0, :) = [];
dyn(dyn.RT < 0.1,:)    = [];


%% Prepare YN data

folder = 'yn_pooled';
load(fullfile(figdir, folder, 'dataPF_joined_and_sep_fits_combined_pffit.mat'), 'pffit')
pffit_yn_pooled = pffit; clear pffit

% Accumulate parameters of subjects in matrices (2 event conditions)
[threshold_matrix, slope_matrix, lapse_matrix] = deal(nan(size(pffit_yn_pooled,2), 4));
[pc_mat] = deal(nan(size(pffit_yn_pooled,2), 3, 8));
for isubj = 1:length(pffit_yn_pooled)
    for icond = 1:length(pffit_yn_pooled{isubj})
        threshold_matrix(isubj, icond) = pffit_yn_pooled{isubj}{icond}.par(1);
        slope_matrix(isubj, icond)  = pffit_yn_pooled{isubj}{icond}.par(2);
        lapse_matrix(isubj, icond)  = pffit_yn_pooled{isubj}{icond}.par(4);
        pc_mat(isubj, icond, :)     = pffit_yn_pooled{isubj}{icond}.perCor;
    end % condition
end % subject

pc_GA_yn_pool = squeeze(mean(pc_mat,1));
pc_SE_yn_pool = squeeze(std(pc_mat)) ./ sqrt(size(pc_mat,1));






%% Preprocess 2IFC data

load(fullfile(data_dir, 'd701to727_2ifc.mat'))
folder = '2ifc';
load(fullfile(figdir, folder, 'dataPF_joined_and_sep_fits_combined_pffit.mat'), 'pffit')
pffit_2ifc = pffit; clear pffit
d2ifc = dall; clear dall

% Accumulate parameters of subjects in matrices
[threshold_matrix, slope_matrix, lapse_matrix] = deal(nan(size(pffit_2ifc,2), 3));
[pc_mat, numcor_mat, outof_mat] = deal(nan(size(pffit_2ifc,2), 3, 8));
for isubj = 1:length(pffit_2ifc)
    for icond = 1:length(pffit_2ifc{isubj})
        threshold_matrix(isubj, icond) = pffit_2ifc{isubj}{icond}.par(1);
        slope_matrix(isubj, icond)  = pffit_2ifc{isubj}{icond}.par(2);
        lapse_matrix(isubj, icond)  = pffit_2ifc{isubj}{icond}.par(4);
        if isubj == 1
            idvect = [1,2,3,4,5,6,8];
        else
            idvect = 1:8;
        end
        pc_mat(isubj, icond, idvect)     = pffit_2ifc{isubj}{icond}.perCor;
    end % condition
end % subject
pc_GA_2ifc       = squeeze(nanmean(pc_mat,1));
pc_SE_2ifc       = squeeze(nanstd(pc_mat)) ./ sqrt(size(pc_mat,1));




%% Define PF

PF = @PAL_Weibull;

soavect = unique(dyn.soa);
soavect(6) = [];
StimLevels = repmat(soavect', [3,1]);

col.indiv  = {[0.6 1 0.6], [0.6 0.6 1], [1 0.6 0.6], [0.6 0.6 0.6]};
col.points = {[0 0.6 0], 'b', 'r', 'k'};
col.ipML   = 'm';
col.ipBA   = [0.6 0.1 0.1];
col.cond   = {[0 0.6 0], 'b', 'r', 'k'};
cols_1F    = [0.3 0.5 0.3];
cols_1F1S  = [0.3 0.3 0.5];
jitterAmount = 0.01;

lw   = 1;  % PF line widths
dtsz = 5;   % size of PF dots




%% Figure YN pooled and 2IFC psychometric functions

fh3 = figure('color', [1 1 1], 'Position', [0, 0, 350*1.5, 350]);
ha = tight_subplot(2, 3,[0.02 0.014],[0.02],[0.02]);


% Plot YN pooled PFs
nj = 20;
ni = 3; 
id = 1:3;
for i = 1:ni % condition
    axes(ha(id(i)));
    
    for j = 1:nj % subj
        soa  = pffit_yn_pooled{j}{i}.soa;
        soal = round((soa*1000))/1000; % round for plotting to two decimal places
        
        % Plot PFs (individual subject)
        soaHR = min(soa):max(soa)/1000:max(soa);
        pffit_yn_pooledML = PF(pffit_yn_pooled{j}{i}.par,    soaHR);
        plot(soaHR, pffit_yn_pooledML, '-', 'Color', col.indiv{i}, 'LineWidth', lw); hold on
    end % end j
    % Plot data
    errorbar(StimLevels(i,:), pc_GA_yn_pool(i,:), pc_SE_yn_pool(i,:), pc_SE_yn_pool(i,:), 'Color', 'k', 'LineWidth', 1.5); hold on
    plot(StimLevels(i,:), pc_GA_yn_pool(i,:), '-', 'Color', col.points{i}, 'linewidth', 1.5); hold on
    plot(StimLevels(i,:), pc_GA_yn_pool(i,:), 'ko', 'MarkerSize', dtsz, 'MarkerFaceColor', col.points{i}); hold on
    
    % plot settings
    set(gca, 'Xtick', soa); yl = [0.5-0.1, 1.05]; ylim(yl);
    set(gca, 'Xlim', [0 max(soa)+min(soa)])
    set(gca, 'YTickLabel', []);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTick', 0.5:0.25:1.5);
    xlim([min(soa)-0.025, max(soa)+0.025]);
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
end % i


% Plot 2IFC psychometric functions
id = 4:6;
for i = 1:ni % condition
    axes(ha(id(i)));
    
    for j = 1:nj % subj
        soa  = pffit_2ifc{j}{i}.soa;
        soal = round((soa*1000))/1000; % round for plotting to two decimal places
        
        % Plot PFs (individual subject)
        soaHR = min(soa):max(soa)/1000:max(soa);
        pffit_2ifcML = PF(pffit_2ifc{j}{i}.par,    soaHR);
        plot(soaHR, pffit_2ifcML, '-', 'Color', col.indiv{i}, 'LineWidth', lw); hold on
    end % end j
    % Plot data
    errorbar(StimLevels(i,:), pc_GA_2ifc(i,:), pc_SE_2ifc(i,:), pc_SE_2ifc(i,:), 'Color', 'k', 'LineWidth', 1.5); hold on
    plot(StimLevels(i,:), pc_GA_2ifc(i,:), '-', 'Color', col.points{i}, 'linewidth', 1.5); hold on
    plot(StimLevels(i,:), pc_GA_2ifc(i,:), 'ko', 'MarkerSize', dtsz, 'MarkerFaceColor', col.points{i}); hold on
    
    % plot settings
    set(gca, 'Xtick', soa); yl = [0.5-0.1, 1.05]; ylim(yl);
    set(gca, 'Xlim', [0 max(soa)+min(soa)])
    set(gca, 'YTickLabel', []);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTick', 0.5:0.25:1.5);
    xlim([min(soa)-0.025, max(soa)+0.025]);
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02, 0.02])
end % i

fh3.Renderer = 'painters'; 
saveas(fh3, 'D:\dfi_experiment_figures\Paper_figures\iAF\iAF_betw_betabinom\pfs_pool_end_2ifc_svg.svg')



% eof















