% Prepares signal detection parameters (no PFs) for yes-no
% and 2IFC data for the analyses and figure script
%
% Parent script(s): 
%   None
%   
% Children script(s): 
%   dfi_sdparams_by_sound_pairing.m
%
% Sibling script(s):
%   None
%
%
% DETAILS
%
%       Stimulus overview
%
%     STIM1   Stim2      ID
%       V1    V1V2   =   2            V 
%      V1V2    V1    =   3            V 
%       A1    A1A2   =   4            A only
%      V1A1   V2A1   =   5            Fus
%      V2A1   V1A1   =   6            Fus 
%      A1A2    A1    =   7            A 
%      V1A2   V2A2   =   8            Fis
%      V2A2   V1A2   =   9            Fis 
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
% Last modified 02/05/2021


restoredefaultpath
clear all
close all


%% *** Get 2ifc data ***

% Paths, and data
addpath(genpath('dfi'));
dfi_startup
data_dir = fullfile('dfi_experiment_data', 'data', 'experiment');
load(fullfile(data_dir, 'd701to727_2ifc.mat'))


%% *** Preprocessing and participant selection ***
% delete no response and bad response trials
draw  = dall;
dall(dall.resp  == 0, :) = [];
dall(dall.badRT ~= 0, :) = [];
dall(dall.RT < 0.1, :)   = [];

% pool over equivalent conditions (eg 1F 2F and 2F 1F)
% 2 and 3 / v only
% 5 and 6 / fusion
% 8 and 9 / fission
dall.trlid(dall.trlid == 2) = 3;
dall.trlid(dall.trlid == 5) = 6;
dall.trlid(dall.trlid == 8) = 9;
dall.trlid(dall.trlid == 4) = 7;

% get descriptive statistics of remaining data
dstats = tabulate(dall.partid);
dstats_distr = dstats(dstats(:,2) > 0,2);
mean(dstats_distr)
max(dstats_distr)
min(dstats_distr)

% get some descriptive statistics on multiple response trials
ddel = draw;
ddel = ddel(ddel.badRT ~= 0, :);
tabulate(ddel.partid)
ddel_stats = tabulate(ddel.partid);
ddel_del_distr = ddel_stats(dstats(:,2) > 0,2);
mean(ddel_del_distr)
max(ddel_del_distr)
min(ddel_del_distr)

% get some descriptive statistics on deleted trials due to fast RTs
dtemp = draw;
dtemp(dtemp.resp  == 0, :) = [];
dtemp(dtemp.badRT ~= 0, :) = [];

dRT = dtemp(dtemp.RT < 0.1,:);

drt_stats = tabulate(dRT.partid);
rt_del_distr = drt_stats(dstats(1:end-1,2) > 0,2);
mean(rt_del_distr)
max(rt_del_distr)
min(rt_del_distr)


%% SIGNAL DETECTION ANALYSES
dall = draw;
dall(dall.resp  == 0, :) = [];
dall(dall.badRT ~= 0, :) = [];
dall(dall.RT < 0.1, :)   = [];

d = dall;

subjvect = unique(dall.partid);
    
clear fh
fh = figure('color', [1 1 1], 'Position', [100, 100, 1800, 900]);
ha = tight_subplot(5,4,[.04 .005],[.1 .05],[.05 .025]);

soaAll = unique(dall.soa);
soavect= soaAll;
nsoa   = length(soaAll);
pC_mat = nan(nsoa, 3, length(subjvect));
dP_mat = nan(nsoa, 3, length(subjvect));
C_mat  = nan(nsoa, 3, length(subjvect));
dP_adj = nan(nsoa, 3, length(subjvect));
C_adj  = nan(nsoa, 3, length(subjvect));
pC_adj = nan(nsoa, 3, length(subjvect));
dP_unbiased = nan(nsoa, 3, length(subjvect));
N_trls = nan(nsoa, 6, length(subjvect));

% loop over participants
for isubj = 1:length(subjvect)
    
    d = dall(dall.partid == subjvect(isubj),:);
    
    if ~isempty(d)
        
        [sdm23, sdm23_adj] = dfi_SDM_2ifc(d,  3, 2, 1);
        [sdm56, sdm56_adj] = dfi_SDM_2ifc(d,  6, 5, 1);
        [sdm89, sdm89_adj] = dfi_SDM_2ifc(d,  9, 8, 1);
        
        soavect = zeros(nsoa,1);
        for isoa = 1:length(soaAll)
            if sum(soaAll(isoa) == unique(d.soa))
                soavect(isoa) = 1;
            end
        end
        soavect = logical(soavect);
        
        dP_mat(soavect, 1, isubj) = sdm23.dP;
        dP_mat(soavect, 2, isubj) = sdm56.dP;
        dP_mat(soavect, 3, isubj) = sdm89.dP;
        
        pC_mat(soavect, 1, isubj) = sdm23.pC;
        pC_mat(soavect, 2, isubj) = sdm56.pC;
        pC_mat(soavect, 3, isubj) = sdm89.pC;
        
        C_mat(soavect, 1, isubj) = sdm23.C;
        C_mat(soavect, 2, isubj) = sdm56.C;
        C_mat(soavect, 3, isubj) = sdm89.C;
        
        dP_adj(soavect, 1, isubj) = sdm23_adj.dP_adj;
        dP_adj(soavect, 2, isubj) = sdm56_adj.dP_adj;
        dP_adj(soavect, 3, isubj) = sdm89_adj.dP_adj;
        
        pC_adj(soavect, 1, isubj) = sdm23_adj.pC_adj;
        pC_adj(soavect, 2, isubj) = sdm56_adj.pC_adj;
        pC_adj(soavect, 3, isubj) = sdm89_adj.pC_adj;
        
        C_adj(soavect, 1, isubj) = sdm23_adj.C_adj;
        C_adj(soavect, 2, isubj) = sdm56_adj.C_adj;
        C_adj(soavect, 3, isubj) = sdm89_adj.C_adj;
        
        soas = unique(dall.soa);
        for isoa = 1:length(soas)
            cond_idx = [2,3,5,6,8,9];
            for icond = 1:length(cond_idx)
                N_trls(isoa, icond, isubj) = ...
                    sum(dall.soa==soas(isoa) & ...
                    dall.partid==subjvect(isubj) & ...
                    dall.trlid==cond_idx(icond));
            end
        end
        
        % plot
        axes(ha(isubj))
        soa= unique(d.soa);
        soa= round((soa*1000))/1000;
        % Plot d'
        plot(soa(isfinite(sdm23.dP)), sdm23.dP(isfinite(sdm23.dP)), 'ks:', 'color', [0 0.75 0], 'LineWidth', 2);  hold on;
        plot(soa(isfinite(sdm56.dP)), sdm56.dP(isfinite(sdm56.dP)), 'bo:', 'LineWidth', 2)
        plot(soa(isfinite(sdm89.dP)), sdm89.dP(isfinite(sdm89.dP)), 'rd:', 'LineWidth', 2)
        % Plot criterion line
        plot(soa(isfinite(sdm23.C)), sdm23.C(isfinite(sdm23.C)), 'k', 'color', [0 0.75 0], 'LineWidth', 1);
        plot(soa(isfinite(sdm56.C)), sdm56.C(isfinite(sdm56.C)), 'k', 'color', [0 0 1], 'LineWidth', 1);
        plot(soa(isfinite(sdm89.C)), sdm89.C(isfinite(sdm89.C)), 'k', 'color', [1 0 0], 'LineWidth', 1);
        set(gca, 'Xlim', [-0.005 0.26])
        set(gca, 'Ylim', [-2.1 5.1])
        set(gca, 'Xtick', soa, 'fontsize', 14);
        if length(subjvect)-isubj > 3
            set(gca, 'XTickLabel', []);
        else
            xlabels = num2cell(soa);
            xlabels = [xlabels(1:2); {[]}; xlabels(4:end)];
            set(gca,'XTickLabel',xlabels)
            xlabel('SOA (s)', 'fontsize', 14);
            rotateXLabels(gca,45)
        end
        title(sprintf('P %d', subjvect(isubj)), 'fontsize', 14);
        if isubj == 2, text(soa(numel(soa)-1), 0, sprintf('       pH > pF\n\n       pH < pF'), ...
                'HorizontalAlignment', 'left', 'VerticalAlignment'  , 'middle'); end
        if ismember(isubj, [1:4:length(subjvect)])
            ylabel('d-prime', 'fontsize', 14);
        else
            set(gca, 'YTickLabel', []);
        end
        gridxy(get(gca,'xtick'),get(gca,'ytick'),'color',[.8 .8 .8],'linewidth',1, 'linestyle', ':')
    end
end % sessions loop

suptitle(sprintf('Sensitivity, all sessions'))

% Show number of trials by flash-beep pairing, when considering
% the intermediate 4 SOAs:
cond_labels = {'1F0B', '2F0B', '1F1B', '2F1B', '1F2B', '2F2B'};
Nmat = squeeze(sum(N_trls(3:6,:,:)));
Nmat_avg = mean(Nmat,2);
Nmat_se = std(Nmat') ./ sqrt(size(Nmat,2));
for icond = 1:6
    fprintf('\nNumber of trials for %s = %f +- %f', ...
        cond_labels{icond}, Nmat_avg(icond), Nmat_se(icond))
end
disp(' ')
        

% save signal detection paramters
mkdir(fullfile(data_dir, '2IFC'))
save(fullfile( data_dir, '2IFC', 'SD_params.mat' ), ...
    'dP_mat', 'C_mat', 'dP_adj', 'C_adj', 'pC_mat', 'pC_adj', 'soaAll', 'N_trls');
close all










% Paths, and data
load(fullfile(data_dir, 'd701to727_yn.mat'))


%% *** Preprocessing and participant selection ***
% delete no response and bad response trials
draw  = dall;
dall(dall.resp  == 0, :) = [];
dall(dall.badRT ~= 0, :) = [];
dall(dall.RT < 0.1,:)    = [];

% get descriptive statistics of remaining data
dstats = tabulate(dall.partid);
dstats_distr = dstats(dstats(:,2) > 0,2);
mean(dstats_distr)
max(dstats_distr)
min(dstats_distr)

% get some descriptive statistics on multiple response trials
ddel = draw;
ddel = ddel(ddel.badRT ~= 0, :);
tabulate(ddel.partid)
ddel_stats = tabulate(ddel.partid);
ddel_del_distr = ddel_stats(dstats(:,2) > 0,2);
mean(ddel_del_distr)
max(ddel_del_distr)
min(ddel_del_distr)

% get some descriptive statistics on deleted trials due to fast RTs
dtemp = draw;
dtemp(dtemp.resp  == 0, :) = [];
dtemp(dtemp.badRT ~= 0, :) = [];

dRT = dtemp(dtemp.RT < 0.1,:);

drt_stats = tabulate(dRT.partid);
rt_del_distr = drt_stats(dstats(1:end-1,2) > 0,2);
mean(rt_del_distr)
max(rt_del_distr)
min(rt_del_distr)


%% SIGNAL DETECTION ANALYSES
dall_old = dall;
dall(dall.soa == (11/120),:) = [];

subjvect = unique(dall.partid);
    
clear fh
fh = figure('color', [1 1 1], 'Position', [100, 100, 1800, 900]);
ha = tight_subplot(5,4,[.04 .005],[.1 .05],[.05 .025]);

soaAll = unique(dall.soa);
nsoa   = length(soaAll);
pC_mat = nan(nsoa, 3, length(subjvect));
dP_mat = nan(nsoa, 3, length(subjvect));
C_mat  = nan(nsoa, 3, length(subjvect));
numPos = nan(nsoa, 3, length(subjvect));
numTot = nan(nsoa, 3, length(subjvect));
N_trls = nan(nsoa, 6, length(subjvect));


% loop over participants
for isubj = 1:length(subjvect)
    
    d = dall(dall.partid == subjvect(isubj),:);
    
    if ~isempty(d)
        [sdm23, sdm23_adj] = dfi_SDM(d,  3, 2, 1);
        [sdm56, sdm56_adj] = dfi_SDM(d,  6, 5, 1);
        [sdm89, sdm89_adj] = dfi_SDM(d,  9, 8, 1);
        
        soavect = zeros(nsoa,1);
        for isoa = 1:length(soaAll)
            if sum(soaAll(isoa) == unique(d.soa))
                soavect(isoa) = 1;
            end
        end
        soavect = logical(soavect);
        
        dP_mat(soavect, 1, isubj) = sdm23.dP;
        dP_mat(soavect, 2, isubj) = sdm56.dP;
        dP_mat(soavect, 3, isubj) = sdm89.dP;
        
        pC_mat(soavect, 1, isubj) = sdm23.pC;
        pC_mat(soavect, 2, isubj) = sdm56.pC;
        pC_mat(soavect, 3, isubj) = sdm89.pC;
        
        C_mat(soavect, 1, isubj) = sdm23.C;
        C_mat(soavect, 2, isubj) = sdm56.C;
        C_mat(soavect, 3, isubj) = sdm89.C;
        
        dP_adj(soavect, 1, isubj) = sdm23_adj.dP_adj;
        dP_adj(soavect, 2, isubj) = sdm56_adj.dP_adj;
        dP_adj(soavect, 3, isubj) = sdm89_adj.dP_adj;
        
        pC_adj(soavect, 1, isubj) = sdm23_adj.pC_adj;
        pC_adj(soavect, 2, isubj) = sdm56_adj.pC_adj;
        pC_adj(soavect, 3, isubj) = sdm89_adj.pC_adj;
        
        C_adj(soavect, 1, isubj) = sdm23_adj.C_adj;
        C_adj(soavect, 2, isubj) = sdm56_adj.C_adj;
        C_adj(soavect, 3, isubj) = sdm89_adj.C_adj;
        
        soas = unique(dall.soa);
        for isoa = 1:length(soas)
            cond_idx = [2,3,5,6,8,9];
            for icond = 1:length(cond_idx)
                N_trls(isoa, icond, isubj) = ...
                    sum(dall.soa==soas(isoa) & ...
                    dall.partid==subjvect(isubj) & ...
                    dall.trlid==cond_idx(icond));
            end
        end
        
        % plot
        axes(ha(isubj))
        soa= unique(d.soa);
        soa= round((soa*1000))/1000;
        % Plot d'
        plot(soa(isfinite(sdm23.dP)), sdm23.dP(isfinite(sdm23.dP)), 'ks:', 'color', [0 0.75 0], 'LineWidth', 2);  hold on;
        plot(soa(isfinite(sdm56.dP)), sdm56.dP(isfinite(sdm56.dP)), 'bo:', 'LineWidth', 2)
        plot(soa(isfinite(sdm89.dP)), sdm89.dP(isfinite(sdm89.dP)), 'rd:', 'LineWidth', 2)
        % Plot criterion line
        plot(soa(isfinite(sdm23.C)), sdm23.C(isfinite(sdm23.C)), 'k', 'color', [0 0.75 0], 'LineWidth', 1);
        plot(soa(isfinite(sdm56.C)), sdm56.C(isfinite(sdm56.C)), 'k', 'color', [0 0 1], 'LineWidth', 1);
        plot(soa(isfinite(sdm89.C)), sdm89.C(isfinite(sdm89.C)), 'k', 'color', [1 0 0], 'LineWidth', 1);
        
        set(gca, 'Xlim', [-0.005 0.26])
        set(gca, 'Ylim', [-4.1 4.1])
        set(gca, 'Xtick', soa, 'fontsize', 14);
        if length(subjvect)-isubj > 3
            set(gca, 'XTickLabel', []);
        else
            xlabels = num2cell(soa);
            xlabels = [xlabels(1:2); {[]}; xlabels(4:end)];
            set(gca,'XTickLabel',xlabels)
            xlabel('SOA (s)', 'fontsize', 14);
            rotateXLabels(gca,45)
        end
        title(sprintf('P %d', subjvect(isubj)), 'fontsize', 14);
        %if isubj == 1, legend('V1     vs V2', 'V1A1 vs V2A1', 'V1A1 vs V1A2', 'Location', 'SouthEast'); end
        if isubj == 2, text(soa(numel(soa)-1), 0, sprintf('       pH > pF\n\n       pH < pF'), 'HorizontalAlignment', 'left', 'VerticalAlignment'  , 'middle'); end
        if ismember(isubj, [1:4:length(subjvect)])
            ylabel('d-prime', 'fontsize', 14);
        else
            set(gca, 'YTickLabel', []);
        end
        gridxy(get(gca,'xtick'),get(gca,'ytick'),'color',[.8 .8 .8],'linewidth',1, 'linestyle', ':')
    end
end % sessions loop

suptitle(sprintf('Sensitivity, all sessions'))

% Show number of trials by flash-beep pairing, when considering
% the intermediate 4 SOAs:
cond_labels = {'1F0B', '2F0B', '1F1B', '2F1B', '1F2B', '2F2B'};
Nmat = squeeze(sum(N_trls(3:6,:,:)));
Nmat_avg = mean(Nmat,2);
Nmat_se = std(Nmat') ./ sqrt(size(Nmat,2));
for icond = 1:6
    fprintf('\nNumber of trials for %s = %f +- %f', ...
        cond_labels{icond}, Nmat_avg(icond), Nmat_se(icond))
end
disp(' ')

% save signal detection paramters
mkdir(fullfile( data_dir, 'yesno'))
save(fullfile( data_dir, 'yesno', 'SD_params.mat'), ...
    'dP_mat', 'C_mat', 'dP_adj', 'C_adj', 'pC_mat', 'pC_adj', 'soaAll', 'N_trls');


% eof

