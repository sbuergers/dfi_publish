
% Assess PF parameter estimate consistency over tasks and conditions
% (2IFC ~ Yes-no; 2F~2F1B, etc...)
%
% sb, sbuergers@gmail.com, Jan 2020
%

% add project directory to path
addpath(genpath('dfi'))

figdir = 'D:\dfi_experiment_figures\PFs\beta_binom_weibull';


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


function fh = plot_detailed_correlation(mat1, mat2, ttl_labels, xlbl, ylbl)
% ttl_labels (cell-array): titles of each subplot
% xlbl (cell-array): x-label for each subplot
% ylbl (cell-array): y-label for each subplot
    fh = figure('color', 'w', 'position', [50 50 900, 250]);
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
        xlabel(xlbl{icond});
        ylabel(ylbl{icond});
        title(ttl_labels{icond});
    end
end





% eof



