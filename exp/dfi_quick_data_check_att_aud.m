function dfi_quick_data_check_att_aud( data, distr, save_figures )
% dfi_quick_data_check2( data, distr, save_figures )
% Gives a quick idea of how the data look like in the SIFI paradigm,
% specify whether yes-no, or 2IFC task was performed ('yesno', '2IFC').
% distr can be for example 'Normal', or 'logist' ...
%
% Optionally, you can save the figures, it will prompt you where to save
% them, save_figures is a boolean. After setting the path the filenames and
% sub-paths are generated automatically.
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



%% Psychometric function

% delete slow response trials
dall(dall.RT > 1.5, :) = [];
dall(dall.RT < 0.1, :) = [];

% PARAMS
dist       = distr; % Normal, Weibull, Gumbel, Quick, logQuick, HypSec, logist
sg.alpha   = 0.003:0.005:0.15;      % threshold (inflection point)
sg.beta    = 10.^[-1:0.1:3];        % slope (log scale)
sg.gamma   = 0:.025:0.8;            % guess-rate
sg.lambda  = 0.01:.005:0.15;        % lapse-rate
fparams    = [1 1 1 1];             % [alpha, beta, gamma, lambda]

fparamsfis = fparams;

header = '2S trials';
fig_dim = [0 0 700 500];
pffit{1} = dfi_fit_pf(dall(dall.trlid == 7,:), dist, sg, fparams, 0, header, 1, fig_dim);
plotspecs
% Add accuracy in 1 S condition
d = dall(dall.trlid == 4,:);
StimLevels = unique(d.soa);
[NumCorr, OutOfNum] = deal(zeros(size(StimLevels)));
for j=1:length(StimLevels)
    accVect     = d.acc(d.soa == StimLevels(j));
    NumCorr(j)  = sum(accVect);
    OutOfNum(j) = length(accVect);
end
propCorr = NumCorr./OutOfNum;
subplot(2,2,1)
plot(StimLevels, propCorr, 'k*', 'MarkerSize', 7, 'MarkerFaceColor', 'k');

% ****************************

% make plots similar to Samaha 2015, proportion correct per soa over both
% 1S and 2S trials:
header = '1S & 2S trials';
dsamaha = dall;
dsamaha.trlid(dsamaha.trlid == 4) = 7;
sgsamaha = sg;
sgsamaha.gamma = 0.5;
sgfparams = [1 1 0 1];
pffit_sam{1} = dfi_fit_pf(dsamaha(dsamaha.trlid == 7,:), dist, sgsamaha, sgfparams, 0, header, 1, fig_dim);
plotspecs

%% save figures?
if save_figures
    saveas(pffit{1}.fh, fullfile(fig_dir, sprintf('PF_2S_%s.emf', distr)))
    saveas(pffit_sam{1}.fh, fullfile(fig_dir, sprintf('PF_samaha_2S_%s.emf', distr)))
    saveas(fh_data_quality, fullfile(fig_dir, sprintf('Data_quality.emf')))
    saveas(fh_RT_overall, fullfile(fig_dir, sprintf('RT_sounds_only_overall.emf')))
    saveas(fh_RT_by_SOA, fullfile(fig_dir, sprintf('RT_sounds_only_by_soa.emf')))
    close all
end




end % fun

% eof


















