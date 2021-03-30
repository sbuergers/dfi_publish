% Compute performance accuracy for different flash-beep pairings
% in yes-no and yes-no threshold task. 
%
% Parent script(s): 
%   None
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
% Last modified 30/03/2021


%% accuracy for different conditions (1f, 1f1s, ... 2f2s) - YNT

load(fullfile('dfi_experiment_data', 'data', 'experiment', 'd701to727_ynt.mat'))
clear d7*

% delete no response and bad response trials
draw_ynt  = dall;
dynt = draw_ynt;
dynt(dynt.resp  == 0, :) = [];
dynt(dynt.RT < 0.1,:)    = [];
dynt(dynt.trlid == 1, :) = [];
dynt(dynt.badRT ~= 0,:)  = [];


% Accumulate parameters of subjects in matrices (2 event conditions)
partvect = unique(dynt.partid);
condvect = unique(dynt.trlid);
soavect = unique(dynt.soa);
[pc_mat, numcor_mat, outof_mat] = deal(nan(length(partvect), length(condvect)));
for isubj = 1:length(partvect)
    for icond = 1:length(condvect)
            acc = dynt.acc(dynt.partid == partvect(isubj) & ...
                dynt.trlid == condvect(icond));
            pc_mat(isubj,icond) = sum(acc) / length(acc);
            numcor_mat(isubj,icond) = sum(acc);
            outof_mat(isubj,icond) = length(acc);
    end % condition
end % subject

% Accumulate parameters of subjects in matrices (1 event conditions)
condition_vect = [2,5];
[per_cor] = nan(20, 2);
for isubj = 1:20
    for icond = 1:2
        dtemp = dynt(dynt.partid == subjvect(isubj) & dynt.trlid == condition_vect(icond),:);
        per_cor(isubj,icond) = sum(dtemp.acc) / length(dtemp);
    end % condition
end % subject

% Get one matrix including data averaged for intermediate soas for 2event
% contexts:
pc = pc_mat;
pc(:, [1, 3]) = per_cor;

pc_GA = squeeze(nanmean(pc,1));
pc_SE = squeeze(nanstd(pc)) ./ sqrt(size(pc,1));

disp('Yes-no threshold')
disp(condvect')
disp(pc_GA)
disp(pc_SE)



%% accuracy for different conditions (1f, 1f1s, ... 2f2s) - YESNO interm. SOAs

% delete no response and bad response trials
load(fullfile('dfi_experiment_data', 'data', 'experiment', 'd701to727_yn.mat'))
clear d7*

draw_yn  = dall;
dyn = draw_yn;
dyn(dyn.resp  == 0,:)  = [];
dyn(dyn.RT < 0.1,:)    = [];
dyn(dyn.trlid == 1,:)  = [];
dyn(dyn.badRT ~= 0,:)  = [];


% Accumulate parameters of subjects in matrices (2 event conditions)
partvect = unique(dyn.partid);
condvect = unique(dyn.trlid);
soavect = unique(dyn.soa);
[pc_mat, numcor_mat, outof_mat] = deal(nan(length(partvect), length(condvect), 8));
for isubj = 1:length(partvect)
    for icond = 1:length(condvect)
        for isoa = 1:8
            acc = dyn.acc(dyn.partid == partvect(isubj) & ...
                dyn.trlid == condvect(icond) & ...
                dyn.soa == soavect(isoa));
            pc_mat(isubj,icond,isoa) = sum(acc) / length(acc);
            numcor_mat(isubj,icond,isoa) = sum(acc);
            outof_mat(isubj,icond,isoa) = length(acc);
        end
    end % condition
end % subject

% Accumulate parameters of subjects in matrices (1 event conditions)
condition_vect = [2,5];
[per_cor] = nan(20, 2);
for isubj = 1:20
    for icond = 1:2
        dtemp = dyn(dyn.partid == subjvect(isubj) & dyn.trlid == condition_vect(icond),:);
        per_cor(isubj,icond) = sum(dtemp.acc) / length(dtemp);
    end % condition
end % subject

% Get one matrix including data averaged for intermediate soas for 2event
% contexts:
pc = nanmean(pc_mat(:,:,3:6),3);
pc(:,[1,3]) = per_cor;

pc_GA = squeeze(nanmean(pc,1));
pc_SE = squeeze(nanstd(pc)) ./ sqrt(size(pc,1));

disp('Yes-no intermediate SOAs')
disp(condvect')
disp(pc_GA)
disp(pc_SE)


% eof

