% Plot substantial evidence lines in Bayes factor figures at 1/3 and 3
% following 
%
% Schönbrodt, Felix D., and Eric-Jan Wagenmakers. “Bayes Factor Design 
% Analysis: Planning for Compelling Evidence.” Psychonomic Bulletin & 
% Review 25, no. 1 (February 1, 2018): 128–42. 
% https://doi.org/10.3758/s13423-017-1230-y.
%
% This is required by nature human behaviour for publication, even though
% thresholds at log10(BF)=-0.5 and log10(BF)=0.5 felt very natural. We are
% still close with log10(1/3) = -0.4771 and log10(3.0) = 0.4771.
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
% Last modified Sept. 2021


% Create a new directory for this auxiliary figure
main_dir = fullfile('dfi_experiment_figures', 'Paper_figures', ...
                    'schoenbrodt_bf_thresholds');
mkdir(main_dir);


% Plot figure with Schoenbrodt et al. (2018) BF thresholds of moderate
% evidence
fh2 = figure('color', [1 1 1], 'Position', [0, 0, 427, 400]);
ha = tight_subplot(6, 4,[0.02 0.02],[0.02],[0.02]);

axes(ha(1));

xl = [ 0  3];
yl = [-1, 1];

tif = linspace(0, 3, 100);
plot(tif, repmat(log10(1/3), size(tif)), 'color', [0.75 0.25 0.75]); hold on
plot(tif, repmat(log10(3.0), size(tif)), 'color', [0.75 0.25 0.75]);

xlim(xl)
ylim(yl)

box off
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.02])
set(gca,'XColor','k','YColor','k')

fh2.Renderer = 'painters';
saveas(fh2, fullfile(main_dir, 'bf_thresholds_schoenbrodt.svg'))

close all


% eof

