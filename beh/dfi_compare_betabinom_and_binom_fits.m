% Assess parameter estimate consistency for old PF fits
% (with logistic binomial model) and new PF fits (with
% beta binomial Weibull model).
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


%% 2IFC

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
save(fullfile(figdir, folder, 'data_betabinom_in_pffit_format_2IFC.mat'), 'pffit_2ifc');


ifc_beta_threshold = threshold_matrix;
ifc_beta_slope = slope_matrix;
ifc_beta_guess = guess_matrix;
ifc_beta_lapse = lapse_matrix;


% Import binomial parameter estimates
load(fullfile('dfi_experiment_data', 'data', 'experiment', 'data_for_ulrik', ...
    'preselected_data', '2ifc_pffits_preselected'))
pffit_all = pffit_2ifc; 

ifc_threshold = nan(length(pffit_all), length(pffit_all{1}));
ifc_slope = nan(length(pffit_all), length(pffit_all{1}));
ifc_guess = nan(length(pffit_all), length(pffit_all{1}));
ifc_lapse = nan(length(pffit_all), length(pffit_all{1}));
for isubj = 1:20
    for icond = 1:3
        if ~isempty(pffit_all{isubj}{icond})
            ifc_threshold(isubj,icond) = pffit_all{isubj}{icond}.par(1);
            ifc_slope(isubj,icond) = pffit_all{isubj}{icond}.par(2);
            ifc_guess(isubj,icond) = pffit_all{isubj}{icond}.par(3);
            ifc_lapse(isubj,icond) = pffit_all{isubj}{icond}.par(4);
        end
    end
end

fh = plot_correlation(ifc_threshold, ifc_beta_threshold, ...
                      'binomial' ,'beta binomial');
                  
fh = plot_correlation(log(ifc_slope), ifc_beta_slope, ...
                      'binomial' ,'beta binomial');
                  
fh = plot_correlation(ifc_guess, ifc_beta_guess, ...
                      'binomial' ,'beta binomial');
                  
fh = plot_correlation(ifc_lapse, ifc_beta_lapse, ...
                      'binomial' ,'beta binomial');                
                  

%% Yes-no

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
save(fullfile(figdir, folder, 'data_betabinom_in_pffit_format_yesno.mat'), 'pffit_yn');


yn_beta_threshold = threshold_matrix;
yn_beta_slope = slope_matrix;
yn_beta_guess = guess_matrix;
yn_beta_lapse = lapse_matrix;
yn_beta_eta = eta_matrix;


% Import binomial parameter estimates
load(fullfile('dfi_experiment_data', 'data', 'experiment', 'data_for_ulrik', ...
    'preselected_data', 'ynpool_pffits_preselected'))
pffit_all = pffit_ynpool;

yn_threshold = nan(length(pffit_all), length(pffit_all{1}));
yn_slope = nan(length(pffit_all), length(pffit_all{1}));
yn_guess = nan(length(pffit_all), length(pffit_all{1}));
yn_lapse = nan(length(pffit_all), length(pffit_all{1}));
for isubj = 1:20
    for icond = 1:3
        if ~isempty(pffit_all{isubj}{icond})
            yn_threshold(isubj,icond) = pffit_all{isubj}{icond}.par(1);
            yn_slope(isubj,icond) = pffit_all{isubj}{icond}.par(2);
            yn_guess(isubj,icond) = pffit_all{isubj}{icond}.par(3);
            yn_lapse(isubj,icond) = pffit_all{isubj}{icond}.par(4);
        end
    end
end

fh = plot_correlation(yn_threshold, yn_beta_threshold, ...
                      'binomial' ,'beta binomial');
                  
fh = plot_correlation(log(yn_slope), yn_beta_slope, ...
                      'binomial' ,'beta binomial');
                  
fh = plot_correlation(yn_guess, yn_beta_guess, ...
                      'binomial' ,'beta binomial');
                  
fh = plot_correlation(yn_lapse, yn_beta_lapse, ...
                      'binomial' ,'beta binomial'); 



%% Nested functions

function fh = plot_correlation(mat1, mat2, xlbl, ylbl)
    fh = figure('color', 'w', 'position', [50 50 900, 250]);
    cond_labels = {'0S', '1S', '2S'};
    for icond = 1:3
        subplot(1,3,icond)
        [spearRho, pval] = corr(mat1(:,icond), mat2(:,icond), 'type', 'Spearman', 'rows', 'complete');
        [~,b1,b0] = regression(mat1(:,icond), mat2(:,icond), 'one');
        plot(mat1(:,icond), mat2(:,icond), 'bo', 'markersize', 10); hold on
        l1 = min([mat1(:,icond); mat2(:,icond)]); 
        l2 = max([mat1(:,icond); mat2(:,icond)]);
        line([l1  l2], [b0+b1*l1 b0+b1*l2], 'color', 'm');
        xlim([l1, l2]); ylim([l1, l2]);
        addtext(sprintf('\nrho = %.2f\np-value = %.5f\nN = %.2f', spearRho, pval, 20));
        xlabel(xlbl);
        ylabel(ylbl);
        title(cond_labels{icond});
    end
end





% eof
































