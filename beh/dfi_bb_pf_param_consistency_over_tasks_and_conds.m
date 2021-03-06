% Assess psychometric function parameter estimate consistency
% over tasks (2ifc, yesno) and conditions (0-2beeps)
%
% Parent script(s): 
%   dfi_beta_binom_combine_joined_and_sep_fits.m
%   
% Children script(s): 
%   None
%
% Sibling script(s):
%   dfi_bb_pf_threshold_consistency.m (figure script)
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
                  
fh = plot_correlation(ifc_beta_slope, yn_beta_slope, ...
                      '2IFC' ,'yesno');
                  
fh = plot_correlation(ifc_beta_guess, yn_beta_guess, ...
                      '2IFC' ,'yesno');
                  
fh = plot_correlation(ifc_beta_lapse, yn_beta_lapse, ...
                      '2IFC' ,'yesno');               
                  



%% Parameter consistency over conditions

%% 2IFC

% Prepare data
%         A     B
% col 1: 0B vs 1B
% col 2: 0B vs 2B
% col 3: 1B vs 2B
[threshold_A, threshold_B] = deal(nan(size(ifc_beta_threshold)));
[slope_A, slope_B] = deal(nan(size(ifc_beta_threshold)));
[guess_A, guess_B] = deal(nan(size(ifc_beta_threshold)));
[lapse_A, lapse_B] = deal(nan(size(ifc_beta_threshold)));

% threshold
threshold_A(:,1) = ifc_beta_threshold(:,1);  % 0B
threshold_A(:,2) = ifc_beta_threshold(:,1);  % 0B
threshold_A(:,3) = ifc_beta_threshold(:,2);  % 1B

threshold_B(:,1) = ifc_beta_threshold(:,2);  % 1B
threshold_B(:,2) = ifc_beta_threshold(:,3);  % 2B
threshold_B(:,3) = ifc_beta_threshold(:,3);  % 2B

% slope
slope_A(:,1) = ifc_beta_slope(:,1);  % 0B
slope_A(:,2) = ifc_beta_slope(:,1);  % 0B
slope_A(:,3) = ifc_beta_slope(:,2);  % 1B

slope_B(:,1) = ifc_beta_slope(:,2);  % 1B
slope_B(:,2) = ifc_beta_slope(:,3);  % 2B
slope_B(:,3) = ifc_beta_slope(:,3);  % 2B

% guess
guess_A(:,1) = ifc_beta_guess(:,1);  % 0B
guess_A(:,2) = ifc_beta_guess(:,1);  % 0B
guess_A(:,3) = ifc_beta_guess(:,2);  % 1B

guess_B(:,1) = ifc_beta_guess(:,2);  % 1B
guess_B(:,2) = ifc_beta_guess(:,3);  % 2B
guess_B(:,3) = ifc_beta_guess(:,3);  % 2B

% lapse
lapse_A(:,1) = ifc_beta_lapse(:,1);  % 0B
lapse_A(:,2) = ifc_beta_lapse(:,1);  % 0B
lapse_A(:,3) = ifc_beta_lapse(:,2);  % 1B

lapse_B(:,1) = ifc_beta_lapse(:,2);  % 1B
lapse_B(:,2) = ifc_beta_lapse(:,3);  % 2B
lapse_B(:,3) = ifc_beta_lapse(:,3);  % 2B

fh = plot_detailed_correlation(threshold_A, threshold_B, ...
                               {'0 vs 1 beep', '0 vs 2 beep', '1 vs 2 beep'}, ...
                               {'0B', '0B', '1B'}, {'1B', '2B', '2B'});
                           
fh = plot_detailed_correlation(slope_A, slope_B, ...
                               {'0 vs 1 beep', '0 vs 2 beep', '1 vs 2 beep'}, ...
                               {'0B', '0B', '1B'}, {'1B', '2B', '2B'});
                           
fh = plot_detailed_correlation(guess_A, guess_B, ...
                               {'0 vs 1 beep', '0 vs 2 beep', '1 vs 2 beep'}, ...
                               {'0B', '0B', '1B'}, {'1B', '2B', '2B'});
                           
fh = plot_detailed_correlation(lapse_A, lapse_B, ...
                               {'0 vs 1 beep', '0 vs 2 beep', '1 vs 2 beep'}, ...
                               {'0B', '0B', '1B'}, {'1B', '2B', '2B'});

                  
%% Yesno

% Prepare data
%         A     B
% col 1: 0B vs 1B
% col 2: 0B vs 2B
% col 3: 1B vs 2B
[threshold_A, threshold_B] = deal(nan(size(yn_beta_threshold)));
[slope_A, slope_B] = deal(nan(size(yn_beta_threshold)));
[guess_A, guess_B] = deal(nan(size(yn_beta_threshold)));
[lapse_A, lapse_B] = deal(nan(size(yn_beta_threshold)));

% threshold
threshold_A(:,1) = yn_beta_threshold(:,1);  % 0B
threshold_A(:,2) = yn_beta_threshold(:,1);  % 0B
threshold_A(:,3) = yn_beta_threshold(:,2);  % 1B

threshold_B(:,1) = yn_beta_threshold(:,2);  % 1B
threshold_B(:,2) = yn_beta_threshold(:,3);  % 2B
threshold_B(:,3) = yn_beta_threshold(:,3);  % 2B

% slope
slope_A(:,1) = yn_beta_slope(:,1);  % 0B
slope_A(:,2) = yn_beta_slope(:,1);  % 0B
slope_A(:,3) = yn_beta_slope(:,2);  % 1B

slope_B(:,1) = yn_beta_slope(:,2);  % 1B
slope_B(:,2) = yn_beta_slope(:,3);  % 2B
slope_B(:,3) = yn_beta_slope(:,3);  % 2B

% guess
guess_A(:,1) = yn_beta_guess(:,1);  % 0B
guess_A(:,2) = yn_beta_guess(:,1);  % 0B
guess_A(:,3) = yn_beta_guess(:,2);  % 1B

guess_B(:,1) = yn_beta_guess(:,2);  % 1B
guess_B(:,2) = yn_beta_guess(:,3);  % 2B
guess_B(:,3) = yn_beta_guess(:,3);  % 2B

% lapse
lapse_A(:,1) = yn_beta_lapse(:,1);  % 0B
lapse_A(:,2) = yn_beta_lapse(:,1);  % 0B
lapse_A(:,3) = yn_beta_lapse(:,2);  % 1B

lapse_B(:,1) = yn_beta_lapse(:,2);  % 1B
lapse_B(:,2) = yn_beta_lapse(:,3);  % 2B
lapse_B(:,3) = yn_beta_lapse(:,3);  % 2B

fh = plot_detailed_correlation(threshold_A, threshold_B, ...
                               {'0 vs 1 beep', '0 vs 2 beep', '1 vs 2 beep'}, ...
                               {'0B', '0B', '1B'}, {'1B', '2B', '2B'});
                           
fh = plot_detailed_correlation(slope_A, slope_B, ...
                               {'0 vs 1 beep', '0 vs 2 beep', '1 vs 2 beep'}, ...
                               {'0B', '0B', '1B'}, {'1B', '2B', '2B'});
                           
fh = plot_detailed_correlation(guess_A, guess_B, ...
                               {'0 vs 1 beep', '0 vs 2 beep', '1 vs 2 beep'}, ...
                               {'0B', '0B', '1B'}, {'1B', '2B', '2B'});
                           
fh = plot_detailed_correlation(lapse_A, lapse_B, ...
                               {'0 vs 1 beep', '0 vs 2 beep', '1 vs 2 beep'}, ...
                               {'0B', '0B', '1B'}, {'1B', '2B', '2B'});


% eof

