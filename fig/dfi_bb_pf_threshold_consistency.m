% Make correlation plots to assess consistency of PF threshold
% between yes-no and 2IFC task for each sound condition.
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


%% Parameter consistency over tasks

fh = plot_correlation(ifc_beta_threshold, yn_beta_threshold, ...
                      '2IFC' ,'yesno');            
           
% eof

