% Substitute poor PF parameter estimates from joined
% fits with estimates from separate fits (see details).
%
% Parent script(s): 
%   dfi_beta_binom_weibull_analysis_2ifc.m
%   dfi_beta_binom_weibull_analysis_sep_fits_2ifc.m
%   dfi_beta_binom_weibull_analysis_ynpool.m
%   dfi_beta_binom_weibull_analysis_sep_fits_ynpool.m
%   
% Children script(s): 
%   dfi_bb_pf_param_consistency_over_tasks_and_conds.m
%   dfi_bb_fig_beh1_all_PFs.m
%   dfi_compare_betabinom_and_binom_fits.m
%
% Sibling script(s):
%   None
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


addpath(genpath('dfi'))


% data import directory
figdir = fullfile('dfi_experiment_figures', 'PFs', 'beta_binom_weibull');


% some useful variables and behavioural data
ttl = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
       '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
fn = fullfile('dfi_experiment_data', 'data', 'experiment', 'd701to727_2ifc.mat');
load(fn)
subjvect = unique(dall.partid);


%% GET DATA: 2IFC

% bad_2ifc_joined_fits --> 701, 719
% bad_2ifc_sep_fits -----> 719, FUS

% 2IFC - joined fits
folder = '2ifc';
load(fullfile(figdir, folder, 'dataPF.mat'), 'dataPF')

% Prepare data for figures
condvect = [3,6,9];
pffit = cell(20,1);
for isubj = 1:numel(subjvect)
    for icond = 1:3
        pffit{isubj}{icond}.soa = dataPF{isubj}.data.soa;
        pffit{isubj}{icond}.cond = condvect(icond);
        pffit{isubj}{icond}.par = dataPF{isubj}.fit.paramsFitted(icond,:);
        pffit{isubj}{icond}.eta = dataPF{isubj}.fit.eta;
        pffit{isubj}{icond}.perCor = dataPF{isubj}.data.perCor(icond,:);
        pffit{isubj}{icond}.OutOfNum = dataPF{isubj}.data.OutOfNum(icond,:);
        pffit{isubj}{icond}.NumCorr = dataPF{isubj}.data.NumCorr(icond,:);
    end
end
d2ifc_joined = dataPF; clear dataPF


% 2IFC - separate fits
folder = fullfile('2ifc', 'sep_fits');
load(fullfile(figdir, folder, 'dataPF.mat'), 'dataPF')

bad_fits_cell = {'701', '719'};

% Prepare data for figures
condvect = [3,6,9];
[~, bad_fits] = ismember(bad_fits_cell, ttl);
for isubj = 1:bad_fits
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
d2ifc_sep = dataPF; clear dataPF


% Save pffit (after combining sep and joined fits) for plotting in
% dfi_bb_fig_beh1_all_PFs.m
folder = '2ifc';
save(fullfile(figdir, folder, 'dataPF_joined_and_sep_fits_combined_pffit.mat'), 'pffit')



% Accumulate parameters of subjects in matrices
[threshold_matrix, slope_matrix, guess_matrix, lapse_matrix, eta_matrix] = ...
    deal(nan(20, 3));
for isubj = 1:20
    for icond = 1:3
        threshold_matrix(isubj, icond) = d2ifc_joined{isubj}.fit.paramsFitted(icond, 1);
        slope_matrix(isubj, icond)  = d2ifc_joined{isubj}.fit.paramsFitted(icond, 2);
        guess_matrix(isubj, icond)  = d2ifc_joined{isubj}.fit.paramsFitted(icond, 3);
        lapse_matrix(isubj, icond)  = d2ifc_joined{isubj}.fit.paramsFitted(icond, 4);
        eta_matrix(isubj, icond)    = d2ifc_joined{isubj}.fit.eta;
    end % condition
end % subject

% Substitute bad joined fits with separate fits
for icond = 1:3
    threshold_matrix(1,icond) = d2ifc_sep{1,icond}.fit.paramsFitted(1);
    threshold_matrix(15,icond) = d2ifc_sep{15,icond}.fit.paramsFitted(1);
    
    slope_matrix(1,icond) = d2ifc_sep{1,icond}.fit.paramsFitted(2);
    slope_matrix(15,icond) = d2ifc_sep{15,icond}.fit.paramsFitted(2);
    
    guess_matrix(1,icond) = d2ifc_sep{1,icond}.fit.paramsFitted(3);
    guess_matrix(15,icond) = d2ifc_sep{15,icond}.fit.paramsFitted(3);
    
    lapse_matrix(1,icond) = d2ifc_sep{1,icond}.fit.paramsFitted(4);
    lapse_matrix(15,icond) = d2ifc_sep{15,icond}.fit.paramsFitted(4);
    
    eta_matrix(1,icond) = d2ifc_sep{1,icond}.fit.eta;
    eta_matrix(15,icond) = d2ifc_sep{15,icond}.fit.eta;
end

% Substitute bad separate fits with NaN
threshold_matrix(15,2) = NaN;
slope_matrix(15,2) = NaN;
guess_matrix(15,2) = NaN;
lapse_matrix(15,2) = NaN;
eta_matrix(15,2) = NaN;


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


% Save
folder = '2ifc';
saveas(fh, fullfile(figdir, folder, 'PF_params_joined_and_sep_fits_combined.emf'))
save(fullfile(figdir, folder, 'dataPF_joined_and_sep_fits_combined.mat'), ...
    'threshold_matrix', 'slope_matrix', 'guess_matrix', 'lapse_matrix', 'eta_matrix')
close all



%% GET DATA: Yes no

% bad_yn_joined_fits --> 701, 712, 714, 725
% bad_yn_sep_fits -----> 712, FIS

% Yes-no - joined fits
folder = 'yn_pooled';
load(fullfile(figdir, folder, 'dataPF.mat'), 'dataPF')

% Prepare data for figures
condvect = [3,6,9];
pffit = cell(20,1);
for isubj = 1:numel(subjvect)
    for icond = 1:3
        pffit{isubj}{icond}.soa = dataPF{isubj}.data.soa;
        pffit{isubj}{icond}.cond = condvect(icond);
        pffit{isubj}{icond}.par = dataPF{isubj}.fit.paramsFitted(icond,:);
        pffit{isubj}{icond}.eta = dataPF{isubj}.fit.eta;
        pffit{isubj}{icond}.perCor = dataPF{isubj}.data.perCor(icond,:);
        pffit{isubj}{icond}.OutOfNum = dataPF{isubj}.data.OutOfNum(icond,:);
        pffit{isubj}{icond}.NumCorr = dataPF{isubj}.data.NumCorr(icond,:);
    end
end
dyn_joined = dataPF; clear dataPF


% Yes-no - separate fits
folder = fullfile('yn_pooled', 'sep_fits');
load(fullfile(figdir, folder, 'dataPF.mat'), 'dataPF')

bad_fits_cell = {'701', '712', '714', '725'};

% Prepare data for figures
condvect = [3,6,9];
[~, bad_fits] = ismember(bad_fits_cell, ttl);
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
dyn_sep = dataPF; clear dataPF


% Save pffit (after combining sep and joined fits) for plotting in
% dfi_bb_fig_beh1_all_PFs.m
folder = 'yn_pooled';
save(fullfile(figdir, folder, 'dataPF_joined_and_sep_fits_combined_pffit.mat'), 'pffit')


% Accumulate parameters of subjects in matrices
[threshold_matrix, slope_matrix, guess_matrix, lapse_matrix, eta_matrix] = ...
    deal(nan(20, 3));
for isubj = 1:20
    for icond = 1:3
        threshold_matrix(isubj, icond) = dyn_joined{isubj}.fit.paramsFitted(icond, 1);
        slope_matrix(isubj, icond)  = dyn_joined{isubj}.fit.paramsFitted(icond, 2);
        guess_matrix(isubj, icond)  = dyn_joined{isubj}.fit.paramsFitted(icond, 3);
        lapse_matrix(isubj, icond)  = dyn_joined{isubj}.fit.paramsFitted(icond, 4);
        eta_matrix(isubj, icond)    = dyn_joined{isubj}.fit.eta;
    end % condition
end % subject

% Substitute bad joined fits with separate fits
for subjids = [1, 9, 10, 18]
    for icond = 1:3
        threshold_matrix(subjids,icond) = dyn_sep{subjids,icond}.fit.paramsFitted(1);
        slope_matrix(subjids,icond) = dyn_sep{subjids,icond}.fit.paramsFitted(2);
        guess_matrix(subjids,icond) = dyn_sep{subjids,icond}.fit.paramsFitted(3);
        lapse_matrix(subjids,icond) = dyn_sep{subjids,icond}.fit.paramsFitted(4);
        eta_matrix(subjids,icond) = dyn_sep{subjids,icond}.fit.eta;
    end
end

% Substitute bad separate fits with NaN
threshold_matrix(9,3) = NaN;
slope_matrix(9,3) = NaN;
guess_matrix(9,3) = NaN;
lapse_matrix(9,3) = NaN;
eta_matrix(9,3) = NaN;


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


% Save
folder = 'yn_pooled';
saveas(fh, fullfile(figdir, folder, 'PF_params_joined_and_sep_fits_combined.emf'))
save(fullfile(figdir, folder, 'dataPF_joined_and_sep_fits_combined.mat'), ...
    'threshold_matrix', 'slope_matrix', 'guess_matrix', 'lapse_matrix', 'eta_matrix')
close all




% eof

