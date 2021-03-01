function dfi_quick_data_check2( data, paradigm, distr, save_figures )
% dfi_quick_data_check2( data, paradigm, distr, save_figures )
% Gives a quick idea of how the data look like in the SIFI paradigm,
% specify whether yes-no, or 2IFC task was performed ('yesno', '2IFC').
% distr can be for example 'Normal', or 'logist' ...
%
% Optionally, you can save the figures, it will prompt you where to save
% them, save_figures is a boolean. After setting the path the filenames and
% sub-paths are generated automatically.
%
%       Stimulus overview 2IFC
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
%         Stimulus overview yesno
%             V A
%         1   0 0
%         2   1 0
%         3   2 0
%         4   0 1
%         5   1 1     % FUSION control
%         6   2 1     % FUSION
%         7   0 2
%         8   1 2     % FISSION
%         9   2 2     % FISSION control


if ~exist('paradigm', 'var')
    paradigm = 'yesno';
end


if ~exist('distr', 'var')
    distr = 'logist';
end


if ~exist('save_figures', 'var')
    save_figures = 0;
end


if save_figures
    [fig_dir] = uigetdir('E:\dfi_experiment_figures', 'Select path to save figures to');
end



%% *** preprocessing ***

dall = data;

% delete no response and bad response trials
draw  = dall;
dall(dall.resp  == 0, :) = [];
dall(dall.badRT ~= 0, :) = [];



%% *** 2IFC ***

% note that I can use 2IFC for yesno as well, then I'm simply plotting
% figures similar to Samaha et al. (2015). For example, I pool over 1 flash
% and two flash trials and simply plot accuracy. If participants have a
% bias towards responding 1 flash in general, then this will average out
% between the two and we can fit the guess rate to 0.5 assuming that there
% is no discrimination ability at the low SOAs. 
%
% However, this does not work well for the fission condition. In the flash
% only condition people are prone to responding 1 flash for both 1 flash
% and two flash trials for low SOAs. In the fission condition the same is
% true. In the fission condition, however, people are prone to responding 1
% flash for the 2f2s condition and 2 flashes for the fission condition.
% Therefore, a bias towards responding 1 flash overall does not exist and
% a guess rate of 0.5 cannot be assumed.

if strcmp(paradigm, '2IFC')
    
    %% Psychometric functions
    % PARAMS
    dist       = distr; % Normal, Weibull, Gumbel, Quick, logQuick, HypSec, logist
    sg.alpha   = 0.008:0.005:0.15;      % threshold (inflection point)
    sg.beta    = 10.^[0.3010:0.005:3];  % slope (log scale)
    sg.gamma   = 0.5;                   % guess-rate
    sg.lambda  = 0.01:.005:0.15;        % lapse-rate
    fparams    = [1 1 0 1];
    
    % pool over equivalent conditions (eg 1F 2F and 2F 1F)
    % 2 and 3 / v only
    % 5 and 6 / fusion
    % 8 and 9 / fission
    dall.trlid(dall.trlid == 2) = 3;
    dall.trlid(dall.trlid == 5) = 6;
    dall.trlid(dall.trlid == 8) = 9;
    dall.trlid(dall.trlid == 4) = 7;

    header = '1F-2F trials';
    pffit{1} = dfi_fit_pf(dall(dall.trlid == 3,:), dist, sg, fparams, 0, header);
    
    header = '1F1S-2F1S trials';
    pffit{2} = dfi_fit_pf(dall(dall.trlid == 6,:), dist, sg, fparams, 0, header);
    
    header = '1F2S-2F2S trials';
    pffit{3} = dfi_fit_pf(dall(dall.trlid == 9,:), dist, sg, fparams, 0, header);
    
    % Assume same lapse rate for all conditions
    estimFit = 0;
    numRep   = 100;
    outp = dfi_fit_pf_group(dall(ismember(dall.trlid, [3,6,9]),:), dist, sg, fparams, estimFit, numRep);
    
         
    %% Response quality
    fh_data_quality = figure('color', 'w');
    drespPie = [sum(data.resp == 0)  / length(data.resp); ...
        sum(data.badRT ~= 0) / length(data.resp); ...
        sum(data.badRT == 0 & data.resp ~= 0) / length(data.resp)];
    pie(drespPie, {'non-responses', 'multiple responses', 'normal responses'});
    title('Data quality');
    
    
    %% Response times
    
    % response times overall and per condition per SOA
    fh_RT_by_SOA = figure('color', 'w', 'Position', [40 40 1500 950]);
    colormap(prism);
    SOAs = unique(dall.soa);
    Nsoa = length(SOAs);
    for isoa = 1:Nsoa
        subplot(ceil(Nsoa/2),2,isoa)
        [counts1, binCenters1] = hist(dall.RT(dall.trlid == 3 & dall.soa == SOAs(isoa)), 10);
        [counts2, binCenters2] = hist(dall.RT(dall.trlid == 6 & dall.soa == SOAs(isoa)), 10);
        [counts3, binCenters3] = hist(dall.RT(dall.trlid == 9 & dall.soa == SOAs(isoa)), 10);
        plot(binCenters1, counts1, 'g-'); hold on;
        plot(binCenters2, counts2, 'b-');
        plot(binCenters3, counts3, 'r-'); xlim([0 2])
        title(sprintf('SOA %g', round(SOAs(isoa)*1000)/1000));
        if isoa > 6, xlabel('Response time (s)'); end;
        if ismember(isoa, [1:2:Nsoa]), ylabel('Number of responses'); end;
        grid on; plotspecs; 
    end
    suptitle(sprintf('Response times\nParticipant %g, session %s', dall.partid(1), num2str(unique(dall.sess))));
    legend('2F','2F1S','1F2S','location','northeast')
    
    
    % response times overall and per condition
    fh_RT_overall = figure('color', 'w', 'Position', [40 40 600 500]);
    colormap(prism);
    [counts1, binCenters1] = hist(dall.RT(dall.trlid == 3), 10);
    [counts2, binCenters2] = hist(dall.RT(dall.trlid == 6), 10);
    [counts3, binCenters3] = hist(dall.RT(dall.trlid == 9), 10);
    plot(binCenters1, counts1, 'g-'); hold on;
    plot(binCenters2, counts2, 'b-');
    plot(binCenters3, counts3, 'r-'); xlim([0 2])
    title(sprintf('Participant %g, session %s', dall.partid(1), num2str(unique(dall.sess))));
    xlabel('Response time (s)'); ylabel('Number of responses');
    suptitle(sprintf('Response times'));
    plotspecs; grid on;
    legend('2F','2F1S','1F2S','location','northeast')
    
    
    %% save figures?
    if save_figures
        saveas(pffit{1}.fh, fullfile(fig_dir, sprintf('PF_2F_%s.emf', distr)))
        saveas(pffit{2}.fh, fullfile(fig_dir, sprintf('PF_Fus_%s.emf', distr)))
        saveas(pffit{3}.fh, fullfile(fig_dir, sprintf('PF_Fis_%s.emf', distr)))
        saveas(outp.fh, fullfile(fig_dir, sprintf('PFs_%s.emf', distr)))
        saveas(fh_data_quality, fullfile(fig_dir, sprintf('Data_quality.emf')))
        close all
    end
    
    
%% *** yesno ***
    
elseif strcmp(paradigm, 'yesno')
    
    %% Psychometric function
    % PARAMS
    dist       = distr; % Normal, Weibull, Gumbel, Quick, logQuick, HypSec, logist
    sg.alpha   = 0.003:0.005:0.15;      % threshold (inflection point)
    sg.beta    = 10.^[-1:0.1:3];        % slope (log scale)
    sg.gamma   = 0:.025:0.8;            % guess-rate
    sg.lambda  = 0.01:.005:0.15;        % lapse-rate
    fparams    = [1 1 1 1];             % [alpha, beta, gamma, lambda]
    
    fparamsfis = fparams;
    
    header = '2S trials';
    pffit{1} = dfi_fit_pf(dall(dall.trlid == 7,:), dist, sg, fparams, 0, header);
    
    % ****************************
    
    % make plots similar to Samaha 2015, proportion correct per soa over both
    % 1S and 2S trials:
    header = '1S & 2S trials';
    dsamaha = dall;
    dsamaha.trlid(dsamaha.trlid == 4) = 7;
    sgsamaha = sg;
    sgsamaha.gamma = 0.5;
    sgfparams = [1 1 0 1];
    pffit_sam{1} = dfi_fit_pf(dsamaha(dsamaha.trlid == 7,:), dist, sgsamaha, sgfparams, 0, header);
    
    %% Response quality
    fh_data_quality = figure('color', 'w');
    drespPie = [sum(data.resp == 0)  / length(data.resp); ...
        sum(data.badRT ~= 0) / length(data.resp); ...
        sum(data.badRT == 0 & data.resp ~= 0) / length(data.resp)];
    pie(drespPie, {'non-responses', 'multiple responses', 'normal responses'});
    title('Data quality');
    
        
    %% Response times
    
    % response times overall and per condition
    fh_RT_overall = figure('color', 'w', 'Position', [40 40 600 500]);
    colormap(prism);
    [counts1, binCenters1] = hist(dall.RT(dall.trlid == 4), 10);
    [counts2, binCenters2] = hist(dall.RT(dall.trlid == 7), 10);
    plot(binCenters1, counts1, 'g-'); hold on;
    plot(binCenters2, counts2, 'b-');
    title(sprintf('Participant %g, session %s', dall.partid(1), num2str(unique(dall.sess))));
    xlabel('Response time (s)'); ylabel('Number of responses');
    suptitle(sprintf('Response times'));
    plotspecs; grid on;
    legend('1S','2S','location','northeast')
    
    
    % response times overall and per condition per SOA
    fh_RT_by_SOA = figure('color', 'w', 'Position', [40 40 1500 950]);
    colormap(prism);
    SOAs = unique(dall.soa);
    Nsoa = length(SOAs);
    for isoa = 1:Nsoa
        subplot(ceil(Nsoa/2),2,isoa)
        [counts1, binCenters1] = hist(dall.RT(dall.trlid == 4 & dall.soa == SOAs(isoa)), 10);
        [counts2, binCenters2] = hist(dall.RT(dall.trlid == 7 & dall.soa == SOAs(isoa)), 10);
        plot(binCenters1, counts1, 'g-'); hold on;
        plot(binCenters2, counts2, 'b-'); xlim([0 2]);
        title(sprintf('SOA %g', round(SOAs(isoa)*1000)/1000));
        if isoa > 6, xlabel('Response time (s)'); end;
        if ismember(isoa, [1:2:Nsoa]), ylabel('Number of responses'); end;
        grid on; plotspecs; 
    end
    suptitle(sprintf('Response times\nParticipant %g, session %s', dall.partid(1), num2str(unique(dall.sess))));
    legend('1S','2S','location','northeast')

    
    %% save figures?
    if save_figures
        saveas(pffit{1}.fh, fullfile(fig_dir, sprintf('PF_2S_%s.emf', distr)))
        saveas(pffit_sam{1}.fh, fullfile(fig_dir, sprintf('PF_samaha_2S_%s.emf', distr)))
        saveas(fh_data_quality, fullfile(fig_dir, sprintf('Data_quality.emf')))
        saveas(fh_RT_overall, fullfile(fig_dir, sprintf('RT_sounds_only_overall.emf')))
        saveas(fh_RT_by_SOA, fullfile(fig_dir, sprintf('RT_sounds_only_by_soa.emf')))
        close all
    end
    
end % if paradigm



end % fun

% eof


















