% Assess consistency of d' and bias over yes-no and yes-no threshold task
% within subjects.
%
% Parent script(s): 
%   dfi_sdt_analyses.m
%   
% Children script(s): 
%   None
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
% Last modified 08/05/2021


restoredefaultpath
clear all
close all

addpath(genpath('dfi'));
dfi_startup

data_dir = fullfile('dfi_experiment_data', 'data', 'experiment');
fig_dir = fullfile('dfi_experiment_figures', 'Paper_figures', 'iAF', 'beh');
mkdir(fig_dir)


%% -------Yes-no intermediate SOAs--------

%% Statistics on d' and C data

load(fullfile( data_dir, 'yesno', 'SD_params.mat' ));

dP_yesno = dP_adj;
C_yesno = C_adj;

dP_yesno = squeeze(nanmean(dP_yesno(3:6,:,:)));
C_yesno = squeeze(nanmean(C_yesno(3:6,:,:)));


%% -------Yes-no threshold-------

load(fullfile('dfi_experiment_data', 'data', 'experiment', 'd701to727_ynt.mat'))
clear d7*

% delete no response and bad response trials
draw_ynt  = dall;
dynt = draw_ynt;
dynt(dynt.resp  == 0, :) = [];
dynt(dynt.RT < 0.1,:)    = [];
dynt(dynt.trlid == 1, :) = [];
dynt(dynt.badRT ~= 0,:)  = [];


%% Yes-no threshold: d', C

dP_mat = nan(20,3);
C_mat = nan(20,3);
pC_mat = nan(20,3);
actualC = nan(20,3);
partvect = unique(dynt.partid);
for isubj = 1:20

    beh = dynt(dynt.partid == partvect(isubj),:);
    % 3.1) Type 1 signal detection parameters
    %S %N
    [~,sdm23] = dfi_SDM(beh,  3, 2, 0);
    [~,sdm56] = dfi_SDM(beh,  6, 5, 0);
    [~,sdm89] = dfi_SDM(beh,  9, 8, 0);
    
    dP_mat(isubj,1) = sdm23.dP_adj;
    C_mat(isubj,1) = sdm23.C_adj;
    pC_mat(isubj,1) = sdm23.pC_adj;
    
    dP_mat(isubj,2) = sdm56.dP_adj;
    C_mat(isubj,2) = sdm56.C_adj;
    pC_mat(isubj,2) = sdm56.pC_adj;
    
    dP_mat(isubj,3) = sdm89.dP_adj;
    C_mat(isubj,3) = sdm89.C_adj;
    pC_mat(isubj,3) = sdm89.pC_adj;
    
    actualC(isubj,1) = sdm23.actualC;
    actualC(isubj,2) = sdm56.actualC;
    actualC(isubj,3) = sdm89.actualC;
    
end

dP_ynt = dP_mat';
C_ynt = C_mat';


%% Correlate d' and bias over experiments

fh1 = plot_correlation(dP_yesno', dP_ynt', 'yesno' ,'yesno threshold');

fh2 = plot_correlation(C_yesno', C_ynt', 'yesno' ,'yesno threshold');


% eof

