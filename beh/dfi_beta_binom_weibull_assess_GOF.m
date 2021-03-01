% Visualize Goodness-of-fit of beta-binomial 
% psychometric function fits
%
% Parent script(s): 
%   dfi_beta_binom_weibull_analysis_2ifc.m
%   dfi_beta_binom_weibull_analysis_sep_fits_2ifc.m
%   dfi_beta_binom_weibull_analysis_ynpool.m
%   dfi_beta_binom_weibull_analysis_sep_fits_ynpool.m
%   
% Children script(s): 
%   dfi_beta_binom_combine_joined_and_sep_fits.m
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

clc
clear all
close all


addpath(genpath('dfi'))

figdir = fullfile('dfi_experiment_figures', 'PFs', 'beta_binom_weibull');




% 2IFC - joined fits
folder = '2ifc';
load(fullfile(figdir, folder, 'dataPF.mat'), 'dataPF')

fprintf('\n')
ttl = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
       '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
   
fh = figure('color', 'w', 'position', [50, 50, 1600, 900]);
for i = 1:length(dataPF)
    
    fprintf('Subject %i\n', i);
    
    boot = dataPF{i}.Bootstrap.Dev; % dataPF{i}.Bootstrap.LL
    fit =  dataPF{i}.fit.Dev; % dataPF{i}.fit.LL
    pDev = sum(boot > fit) / length(boot);
    
    subplot(4,5,i);

    disp('Goodness of fit:');
    fprintf('Deviance: %6.4f\n', fit);
    fprintf('p-value:  %6.4f\n', pDev);
    % plot D* versus sample data Deviance
    hist(boot, 1000); hold on
    line([fit fit], [0 max(ylim)], 'color', 'r');
    text(fit, max(ylim)-max(ylim)/10, sprintf('Deviance: %g', round(fit*1000)/1000), 'FontSize', 10, 'Color', 'r', ...
                                'HorizontalAlignment', 'left', 'VerticalAlignment' , 'middle');
    text(fit+0.5, max(ylim)-max(ylim)/1.1, sprintf('p-value: %g', round(pDev*1000)/1000), 'FontSize', 10, 'Color', 'r', ...
                                'HorizontalAlignment', 'left', 'VerticalAlignment' , 'middle');   
    title(ttl{i});
end
saveas(fh, fullfile(figdir, sprintf('GOF_joined_%s.emf', folder)))
close all




% Yes-no - joined fits
folder = 'yn_pooled';
load(fullfile(figdir, folder, 'dataPF.mat'), 'dataPF')

fprintf('\n')
ttl = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
       '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
   
fh = figure('color', 'w', 'position', [50, 50, 1600, 900]);
for i = 1:length(dataPF)
    
    fprintf('Subject %i\n', i);
    
    boot = dataPF{i}.Bootstrap.Dev; % dataPF{i}.Bootstrap.LL
    fit =  dataPF{i}.fit.Dev; % dataPF{i}.fit.LL
    pDev = sum(boot > fit) / length(boot);
    
    subplot(4,5,i);

    disp('Goodness of fit:');
    fprintf('Deviance: %6.4f\n', fit);
    fprintf('p-value:  %6.4f\n', pDev);
    % plot D* versus sample data Deviance
    hist(boot, 1000); hold on
    line([fit fit], [0 max(ylim)], 'color', 'r');
    text(fit, max(ylim)-max(ylim)/10, sprintf('Deviance: %g', round(fit*1000)/1000), 'FontSize', 10, 'Color', 'r', ...
                                'HorizontalAlignment', 'left', 'VerticalAlignment' , 'middle');
    text(fit+0.5, max(ylim)-max(ylim)/1.1, sprintf('p-value: %g', round(pDev*1000)/1000), 'FontSize', 10, 'Color', 'r', ...
                                'HorizontalAlignment', 'left', 'VerticalAlignment' , 'middle');   
    title(ttl{i});
end
saveas(fh, fullfile(figdir, sprintf('GOF_joined_%s.emf', folder)))
close all




% Yes-no - separate fits
folder = fullfile('yn_pooled', 'sep_fits');
load(fullfile(figdir, folder, 'dataPF.mat'), 'dataPF')

fprintf('\n')
ttl = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
       '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
   
for icond = 1:3
   fh = figure('color', 'w', 'position', [50, 50, 1600, 900]);
   for i = 1:length(dataPF)

       if ismember(i, [1, 9, 10, 18])

           fprintf('Subject %i\n', i);

           boot = dataPF{i,icond}.Bootstrap.Dev; % dataPF{i}.Bootstrap.LL
           fit =  dataPF{i,icond}.fit.Dev; % dataPF{i}.fit.LL
           pDev = sum(boot > fit) / length(boot);

           subplot(4,5,i);

           disp('Goodness of fit:');
           fprintf('Deviance: %6.4f\n', fit);
           fprintf('p-value:  %6.4f\n', pDev);
           % plot D* versus sample data Deviance
           hist(boot, 1000); hold on
           line([fit fit], [0 max(ylim)], 'color', 'r');
           text(fit, max(ylim)-max(ylim)/10, sprintf('Deviance: %g', round(fit*1000)/1000), 'FontSize', 10, 'Color', 'r', ...
               'HorizontalAlignment', 'left', 'VerticalAlignment' , 'middle');
           text(fit+0.5, max(ylim)-max(ylim)/1.1, sprintf('p-value: %g', round(pDev*1000)/1000), 'FontSize', 10, 'Color', 'r', ...
               'HorizontalAlignment', 'left', 'VerticalAlignment' , 'middle');
           title(ttl{i});
       end
   end
   mkdir(fullfile(figdir, folder))
   saveas(fh, fullfile(figdir, folder, sprintf('GOF_joined_cond%i.emf', icond)))
   close all
end




% 2IFC - separate fits
folder = fullfile('2ifc', 'sep_fits');
load(fullfile(figdir, folder, 'dataPF.mat'), 'dataPF')

fprintf('\n')
ttl = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
       '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};
   
for icond = 1:3
   fh = figure('color', 'w', 'position', [50, 50, 1600, 900]);
   for i = 1:length(dataPF)

       if ismember(i, [1, 15])

           fprintf('Subject %i\n', i);

           boot = dataPF{i,icond}.Bootstrap.Dev; % dataPF{i}.Bootstrap.LL
           fit =  dataPF{i,icond}.fit.Dev; % dataPF{i}.fit.LL
           pDev = sum(boot > fit) / length(boot);

           subplot(4,5,i);

           disp('Goodness of fit:');
           fprintf('Deviance: %6.4f\n', fit);
           fprintf('p-value:  %6.4f\n', pDev);
           % plot D* versus sample data Deviance
           hist(boot, 1000); hold on
           line([fit fit], [0 max(ylim)], 'color', 'r');
           text(fit, max(ylim)-max(ylim)/10, sprintf('Deviance: %g', round(fit*1000)/1000), 'FontSize', 10, 'Color', 'r', ...
               'HorizontalAlignment', 'left', 'VerticalAlignment' , 'middle');
           text(fit+0.5, max(ylim)-max(ylim)/1.1, sprintf('p-value: %g', round(pDev*1000)/1000), 'FontSize', 10, 'Color', 'r', ...
               'HorizontalAlignment', 'left', 'VerticalAlignment' , 'middle');
           title(ttl{i});
       end
   end
   mkdir(fullfile(figdir, folder))
   saveas(fh, fullfile(figdir, folder, sprintf('GOF_joined_cond%i.emf', icond)))
   close all
end


% eof

