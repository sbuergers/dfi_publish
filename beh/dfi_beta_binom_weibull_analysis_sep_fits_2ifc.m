% Fit psychometric functions to 2IFC data
%
% Dependencies:
%   Palamedes toolbox (1.8.0)
%   fitPSFmodel.m
%
% Parent script(s): 
%   None
%   
% Children script(s): 
%   dfi_beta_binom_weibull_assess_GOF.m
%   dfi_beta_binom_combine_joined_and_sep_fits.m
%
% Sibling script(s):
%   dfi_beta_binom_weibull_analysis_2ifc.m
%   dfi_beta_binom_weibull_analysis_ynpool.m
%   dfi_beta_binom_weibull_analysis_sep_fits_ynpool.m
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


clc
clear all
close all


%% *** Get data ***


% require palamedes and add relevant paths
addpath(genpath('dfi'))

% set up general figure directory
figdir = fullfile('dfi_experiment_figures', 'PFs', 'beta_binom_weibull');

% set up folder name of type of analysis
folder = fullfile('2ifc', 'sep_fits');
mkdir(fullfile(figdir, folder));

% Manually set data folder and file name
fn = fullfile('dfi_experiment_data', 'data', 'experiment', 'd701to727_2ifc');
load(fn)

S2_conditions = [3,6,8,9];


%% *** Preprocessing and participant selection ***
dall.trlid(dall.trlid == 2) = 3;
dall.trlid(dall.trlid == 5) = 6;
dall.trlid(dall.trlid == 8) = 9;
dall.trlid(dall.trlid == 4) = 7;


% delete no response and bad response trials
draw  = dall;
dall(dall.trlid == 1, :) = [];
dall(dall.resp  == 0, :) = []; % delete non-response trials
dall(dall.badRT ~= 0, :) = []; % delete multiple response trials
dall(dall.RT < 0.1, :)   = [];



% Fit beta binomial model for all participants
ttl = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
       '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
bad_fits_cell = {'701', '719'};
subjvect = unique(dall.partid);
condvect = [3,6,9];
dataPF = cell(numel(subjvect),numel(condvect));
betaBool = true;
for icond = 1:numel(condvect)
    fprintf('\nCondition %i', condvect(icond))
    for isubj = 1:numel(subjvect)
        if ismember(ttl{isubj}, bad_fits_cell)
            fprintf('\nFitting beta-binomial model for participant %i', subjvect(isubj))
            d = dall(dall.partid == subjvect(isubj) & dall.trlid == condvect(icond),:);
            dataPF{isubj,icond} = fitPSFmodel_sep_fits(d,betaBool);
        end
    end
end
fprintf('\n')
% save PF fits (in dataPF format - from Dave's code)
mkdir(fullfile(figdir, folder))
save(fullfile(figdir, folder, 'dataPF.mat'), 'dataPF')


% Prepare data for figures
condvect = [3,6,9];
[~, bad_fits] = ismember(bad_fits_cell, ttl);
pffit = cell(numel(bad_fits),1);
for isubj = 1:numel(bad_fits)
    for icond = 1:3
        pffit{isubj}{icond}.soa = dataPF{bad_fits(isubj), icond}.data.soa;
        pffit{isubj}{icond}.cond = condvect(icond);
        pffit{isubj}{icond}.par = dataPF{bad_fits(isubj), icond}.fit.paramsFitted;
        pffit{isubj}{icond}.eta = dataPF{bad_fits(isubj), icond}.fit.eta;
        pffit{isubj}{icond}.perCor = dataPF{bad_fits(isubj), icond}.data.perCor;
        pffit{isubj}{icond}.OutOfNum = dataPF{bad_fits(isubj), icond}.data.OutOfNum;
        pffit{isubj}{icond}.NumCorr = dataPF{bad_fits(isubj), icond}.data.NumCorr;
    end
end


% Make Figures
ylab   = {'2F', 'FUS', 'FIS'};
supttl = 'Temporal binding windows (2IFC)'; 
fh     = fig_pf_stacked_nice_figures_no_bay(pffit, ...
                                            ylab, bad_fits_cell, ...
                                            'Weibull', [0 0 1000 800]);
% save figures
mkdir(fullfile(figdir, folder));
saveas(fh, fullfile(figdir, folder, sprintf('PFs.emf')))
close all


% Accumulate parameters of subjects in matrices
[threshold_matrix, slope_matrix, guess_matrix, lapse_matrix, eta_matrix] = ...
    deal(nan(size(pffit,2), 3));
[pc_mat, numcor_mat, outof_mat] = deal(nan(size(pffit,2), 3, 7));
for isubj = 1:length(pffit)
    for icond = 1:length(pffit{isubj})
        threshold_matrix(isubj, icond) = pffit{isubj}{icond}.par(1);
        slope_matrix(isubj, icond)  = pffit{isubj}{icond}.par(2);
        guess_matrix(isubj, icond)  = pffit{isubj}{icond}.par(3);
        lapse_matrix(isubj, icond)  = pffit{isubj}{icond}.par(4);
        eta_matrix(isubj, icond)    = pffit{isubj}{icond}.eta;
        pc_mat(isubj, icond, :)     = pffit{isubj}{icond}.perCor(1:7);
        numcor_mat(isubj, icond, :) = pffit{isubj}{icond}.NumCorr(1:7);
        outof_mat(isubj, icond, :)  = pffit{isubj}{icond}.OutOfNum(1:7);
    end % condition
end % subject
pc_GA       = squeeze(mean(pc_mat,1));
pc_SE       = squeeze(std(pc_mat)) ./ sqrt(size(pc_mat,1));
NumCorr_GA  = squeeze(mean(numcor_mat,1));
OutOfNum_GA = squeeze(mean(outof_mat,1));




% Plot parameter distributions
fh = figure('color', 'w', 'Position', [30 30 1200 800]);
x = 1:size(threshold_matrix,1);

% threshold
subplot(231); 
h1 = scatter(x-0.15, threshold_matrix(:,1), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 1 0]); view(90,90); xlim([0 length(threshold_matrix(:,1))+1]); 
set(gca, 'XTick', 1:1:length(threshold_matrix(:,1))); set(gca, 'YTick', 0.008:0.05:0.25); hold on 
h2 = scatter(x, threshold_matrix(:,2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0 1]); view(90,90); xlim([0 length(threshold_matrix(:,1))+1]);
h3 = scatter(x+0.15, threshold_matrix(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0]); view(90,90); xlim([0 length(threshold_matrix(:,1))+1]);
gridxy(get(gca,'xtick'),get(gca,'ytick'),'color',[.8 .8 .8],'linewidth',1, 'linestyle', '-'); ylim([0, 0.25])
xlabel('Participant id'); ylabel('SOA (ms)'); title('Threshold'); set(gca, 'XTickLabel', subjvect); 

% slope
subplot(232); 
h1 = scatter(x-0.15, slope_matrix(:,1), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 1 0]); view(90,90); xlim([0 length(slope_matrix(:,1))+1]); 
set(gca, 'XTick', 1:1:length(slope_matrix(:,1))); hold on 
h2 = scatter(x, slope_matrix(:,2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0 1]); view(90,90); xlim([0 length(slope_matrix(:,1))+1]);
h3 = scatter(x+0.15, slope_matrix(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0]); view(90,90); xlim([0 length(slope_matrix(:,1))+1]);
gridxy(get(gca,'xtick'),get(gca,'ytick'),'color',[.8 .8 .8],'linewidth',1, 'linestyle', '-')
xlabel('Participant id'); ylabel('a.u.'); title('Slope'); set(gca, 'XTickLabel', subjvect); 

% guess rate
subplot(234); 
h1 = scatter(x-0.15, guess_matrix(:,1), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 1 0]); view(90,90); xlim([0 length(guess_matrix(:,1))+1]); 
set(gca, 'XTick', 1:1:length(guess_matrix(:,1))); set(gca, 'YTick', 0:0.05:0.5); hold on 
h2 = scatter(x, guess_matrix(:,2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0 1]); view(90,90); xlim([0 length(guess_matrix(:,1))+1]);
h3 = scatter(x+0.15, guess_matrix(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0]); view(90,90); xlim([0 length(guess_matrix(:,1))+1]);
gridxy(get(gca,'xtick'),get(gca,'ytick'),'color',[.8 .8 .8],'linewidth',1, 'linestyle', '-')
xlabel('Participant id'); ylabel('guess rate (%)'); title('guess rate'); set(gca, 'XTickLabel', subjvect); 
set(findall(gcf,'-property','FontSize'),'FontSize',12)

% lapse rate
subplot(235); 
h1 = scatter(x-0.15, lapse_matrix(:,1), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 1 0]); view(90,90); xlim([0 length(lapse_matrix(:,1))+1]); 
set(gca, 'XTick', 1:1:length(lapse_matrix(:,1))); set(gca, 'YTick', 0:0.025:0.15); hold on 
h2 = scatter(x, lapse_matrix(:,2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0 1]); view(90,90); xlim([0 length(lapse_matrix(:,1))+1]);
h3 = scatter(x+0.15, lapse_matrix(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0]); view(90,90); xlim([0 length(lapse_matrix(:,1))+1]);
gridxy(get(gca,'xtick'),get(gca,'ytick'),'color',[.8 .8 .8],'linewidth',1, 'linestyle', '-')
xlabel('Participant id'); ylabel('Lapse rate (%)'); title('Lapse rate'); set(gca, 'XTickLabel', subjvect); 
set(findall(gcf,'-property','FontSize'),'FontSize',12)

% eta
subplot(236); 
h1 = scatter(x-0.15, eta_matrix(:,1), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 1 0]); view(90,90); xlim([0 length(eta_matrix(:,1))+1]); 
set(gca, 'XTick', 1:1:length(eta_matrix(:,1))); set(gca, 'YTick', 0:0.025:0.15); hold on 
h2 = scatter(x, eta_matrix(:,2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0 1]); view(90,90); xlim([0 length(eta_matrix(:,1))+1]);
h3 = scatter(x+0.15, eta_matrix(:,3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0]); view(90,90); xlim([0 length(eta_matrix(:,1))+1]);
gridxy(get(gca,'xtick'),get(gca,'ytick'),'color',[.8 .8 .8],'linewidth',1, 'linestyle', '-')
xlabel('Participant id'); ylabel('eta rate'); title('eta rate'); set(gca, 'XTickLabel', subjvect); 
set(findall(gcf,'-property','FontSize'),'FontSize',12)


% Save bootstrapped data
saveas(fh, fullfile(figdir, folder, 'PF_params.emf'))


% eof
