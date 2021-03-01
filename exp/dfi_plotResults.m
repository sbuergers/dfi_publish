function [ r ] = dfi_plotResults( s, d )
%dfi_plotResults
%
%INPUT:
% data 'd' from a double flash illusion/fusion paradigm and info-structure 
% 's'.
%
%OUTPUT:
% graphical depiction of response time and accuracy d. Also adds some
% descriptive statistics to the dataset and gives two separate
% datasets for fission and fusion (d12 and d21). 
%
% TRIAL MARKER OVERVIEW
% Trlid:    VA
%           00  -  1
%           10  -  2
%           20  -  3
%           01  -  4
%           11  -  5
%           21  -  6
%           02  -  7
%           12  -  8
%           22  -  9
%--------------------------------------------------------------------------
% SB, July 2015
% sbuergers@gmail.com
%

% Check what type of trials you want to analyze
trltypes = unique(d.trlid);
if isempty(setdiff(trltypes, [1 2 3]))
    fprintf(['Plotting results for Flash stimuli only dataset\n', ...
             '--> Conditions: V1A0, V2A0 as well as A0V0\n ...\n'])
elseif isempty(setdiff(trltypes, [2 3]))
    fprintf(['Plotting results for Flash stimuli only dataset\n', ...
             '--> Conditions: V1A0, V2A0; note there is no V0A0...'])
elseif isempty(setdiff(trltypes, [1 4 7]))
    fprintf(['Plotting results for Sound stimuli only dataset\n', ...
             '--> Conditions: V0A1, V0A2 as well as A0V0\n ...\n'])
elseif isempty(setdiff(trltypes, [4 7]))
    fprintf(['Plotting results for Sound stimuli only dataset\n', ...
             '--> Conditions: V0A1, V0A2, note there is no A0V0\n ...\n'])
elseif isempty(setdiff(trltypes, [1 5 6 8 9]))
    fprintf(['Plotting results for a multimodal fission/fusion paradigm\n', ...
             '--> Conditions: V1A1, V1A2, V2A1, V2A2, as well as A0V0\n ...\n'])
elseif isempty(setdiff(trltypes, [5 6 8 9]))
    fprintf(['Plotting results for a multimodal fission/fusion paradigm\n', ...
             '--> Conditions: V1A1, V1A2, V2A1, V2A2, note there is no A0V0\n ...\n'])
else
    warning('dfi_plotResults:UnexpectedTrialIDs', ...
            'The script is not designed to deal with your configuration of trial types!')
end

% All of these accuracies pertain to key configuration 'AB'
if strcmp(s.paradigm, 'yesno')
    % add accuracy column to data
    acc = zeros(size(d(:,1)));
    acc(d.resp == 1 & ismember(d.trlid, [2,  5,  8,    4])) = 1;
    acc(d.resp == 2 & ismember(d.trlid, [3,  6,  9,    7])) = 1;
    d.acc = acc;
    if strcmp(s.cond, 'multiAttAud')
        acc = zeros(size(d(:,1)));
        acc(d.resp == 1 & ismember(d.trlid, [4 5 6])) = 1;
        acc(d.resp == 2 & ismember(d.trlid, [7 8 9])) = 1;
        d.acc = acc;
    end
elseif strcmp(s.paradigm, '2AFC')
    acc = zeros(size(d(:,1)));
    acc(d.resp == 2 & ismember(d.trlid, [2,  5,  8,    4])) = 1;
    acc(d.resp == 1 & ismember(d.trlid, [3,  6,  9,    7])) = 1;
    d.acc = acc;
elseif strcmp(s.paradigm, 'YN_threshold')
    acc = zeros(size(d(:,1)));
    acc(d.resp < 5 & d.resp > 0 & ismember(d.trlid, [2,  5,  8,    4])) = 1;
    acc(d.resp > 4 &              ismember(d.trlid, [3,  6,  9,    7])) = 1;
    d.acc = acc;
    d.conf = d.resp;
    d.conf(d.resp > 4) = abs(5 - (d.conf(d.resp > 4) - 4));
end

% When key configuration is 'BA' accuracy is simply flipped
d.acc(strcmp(d.rkeycfg, 'BA')) = ~d.acc(strcmp(d.rkeycfg, 'BA'));

% make separate datasets for fusion and fission
d21 = d(d.trlid == 6, :);
d12 = d(d.trlid == 8, :);
% and their controls
d11 = d(d.trlid == 5, :);
d22 = d(d.trlid == 9, :);

% do the same for unimodal conditions
d10 = d(d.trlid == 2, :); % V1
d20 = d(d.trlid == 3, :); % V2
d01 = d(d.trlid == 4, :); % A1
d02 = d(d.trlid == 7, :); % A2


% calculate some simple statistics (overall)
stats.trlid = grpstats(d,'trlid',{'mean', 'std'},'DataVars',{'RT', 'soa', 'acc'});
stats.soa   = grpstats(d,'soa',  {'mean', 'std'},'DataVars',{'RT', 'acc'});
disp('Group statistics (all):')
disp('Grouping Variable: trlid')
disp(stats.trlid)
disp('Grouping Variable: soa')
disp(stats.soa)


%% do the same separately for fission and fusion ...
if any(ismember(d.trlid, [5 6 8 9])) % multimodal condition
    % fission
    stats12.trlid = grpstats(d12,'trlid',{'mean', 'std'},'DataVars',{'RT', 'soa', 'acc'});
    stats12.soa   = grpstats(d12,'soa',  {'mean', 'std'},'DataVars',{'RT', 'acc'});
    disp('Group statistics (FISSION - V1A2 only):')
    disp('Grouping Variable (FISSION - V1A2 only): trlid')
    disp(stats12.trlid)
    disp('Grouping Variable (FISSION - V1A2 only): soa')
    disp(stats12.soa)
    % fission control
    stats22.trlid = grpstats(d22,'trlid',{'mean', 'std'},'DataVars',{'RT', 'soa', 'acc'});
    stats22.soa   = grpstats(d22,'soa',  {'mean', 'std'},'DataVars',{'RT', 'acc'});
    disp('Group statistics (FISSION Control - V2A2 only):')
    disp('Grouping Variable (FISSION Control - V2A2 only): trlid')
    disp(stats22.trlid)
    disp('Grouping Variable (FISSION Control - V2A2 only): soa')
    disp(stats22.soa)

    % fusion
    stats21.trlid = grpstats(d21,'trlid',{'mean', 'std'},'DataVars',{'RT', 'soa', 'acc'});
    stats21.soa   = grpstats(d21,'soa',  {'mean', 'std'},'DataVars',{'RT', 'acc'});
    disp('Group statistics (FUSION - V2A1 only):')
    disp('Grouping Variable: trlid')
    disp(stats21.trlid)
    disp('Grouping Variable: soa')
    disp(stats21.soa)
    % fusion control
    stats11.trlid = grpstats(d11,'trlid',{'mean', 'std'},'DataVars',{'RT', 'soa', 'acc'});
    stats11.soa   = grpstats(d11,'soa',  {'mean', 'std'},'DataVars',{'RT', 'acc'});
    disp('Group statistics (FUSION Control - V1A1 only):')
    disp('Grouping Variable (FUSION Control - V1A1 only): trlid')
    disp(stats11.trlid)
    disp('Grouping Variable (FUSION Control - V1A1 only): soa')
    disp(stats11.soa)
    
    % PLOTS
    disp('Plotting results...')

    % RT (boxplots)
    rtplot = figure;
    subplot(2,3,1)                          % all
    boxplot(d.RT,d.soa)
    title('Response times (both FISSION and FUSION)')
    xlabel('soa'); ylabel('RT');
    subplot(2,3,2)                          % fusion 
    boxplot(d21.RT,d21.soa)
    title('Response times (FUSION: only V2A1)')
    xlabel('soa'); ylabel('RT');
    subplot(2,3,3)                          % fission
    boxplot(d12.RT,d12.soa)
    title('Response times (FISSION: only V1A2)')
    xlabel('soa'); ylabel('RT');
    subplot(2,3,5)                          % fusion control
    boxplot(d11.RT,d11.soa)
    title('Response times (FUSION Control: only V1A1)')
    xlabel('soa'); ylabel('RT');
    subplot(2,3,6)                          % fission control
    boxplot(d22.RT,d22.soa)
    title('Response times (FISSION Control: only V2A2)')
    xlabel('soa'); ylabel('RT');
    set(gcf,'color','w');
    % set(find(0, 'Type', 'title'), 'FontSize', 18)


    % ACC 
    % 2 sd error bars - assuming normpdf, just to get an idea
    accplot = figure;
    subplot(2,3,1)                          % all
    scatter(stats.soa.soa, stats.soa.mean_acc)
    errorbar(stats.soa.soa, stats.soa.mean_acc, stats.soa.std_acc)
    title('Accuracy (both FISSION and FUSION)')
    xlabel('soa'); ylabel('acc')
    subplot(2,3,2)                          % fusion
    scatter(stats21.soa.soa, stats21.soa.mean_acc)
    errorbar(stats21.soa.soa, stats21.soa.mean_acc, stats21.soa.std_acc)
    title('Accuracy (FUSION: only V2A1)')
    xlabel('soa'); ylabel('acc')
    subplot(2,3,3)                          % fission
    scatter(stats12.soa.soa, stats12.soa.mean_acc)
    errorbar(stats12.soa.soa, stats12.soa.mean_acc, stats12.soa.std_acc)
    title('Accuracy (FISSION: only V1A2)')
    xlabel('soa'); ylabel('acc')
    subplot(2,3,5)                          % fusion control
    scatter(stats11.soa.soa, stats11.soa.mean_acc)
    errorbar(stats11.soa.soa, stats11.soa.mean_acc, stats11.soa.std_acc)
    title('Accuracy (FUSION Control: only V1A1)')
    xlabel('soa'); ylabel('acc')
    subplot(2,3,6)                          % fission control
    scatter(stats22.soa.soa, stats22.soa.mean_acc)
    errorbar(stats22.soa.soa, stats22.soa.mean_acc, stats22.soa.std_acc)
    title('Accuracy (FISSION Control: only V2A2)')
    xlabel('soa'); ylabel('acc')
    set(gcf,'color','w'); 
    try
        % Signal detection parameters
        fprintf('\n\nCalculating sensitivity and bias (MULTISENSORY):\n')
        fprintf('\nSignal=V1, Noise=V2:\n');
        sdm56 = dfi_SDM(d, 5, 6); disp(sdm56); % V1A1 - signal, V2A1 - noise
        fprintf('\nSignal=V2, Noise=V1:\n');
        sdm98 = dfi_SDM(d, 9, 8); disp(sdm98); % V2A2 - signal, V1A2 - noise
    catch
        fprintf('\n\nProblem with calculating signal detection parameters. Check Palamedes installation.\n');
    end

end % end multimodal


%% ... or for the unimodal condition at hand
if any(ismember(d.trlid, [2 3 4 7])) % unimodal conditions
    % V1
    stats10.trlid = grpstats(d10,'trlid',{'mean', 'std'},'DataVars',{'RT', 'soa', 'acc'});
    stats10.soa   = grpstats(d10,'soa',  {'mean', 'std'},'DataVars',{'RT', 'acc'});
    disp('Grouping Variable (unimodal - V1 only): trlid')
    disp(stats10.trlid)
    disp('Grouping Variable (unimodal - V1 only): soa')
    disp(stats10.soa)
    % V2
    stats20.trlid = grpstats(d20,'trlid',{'mean', 'std'},'DataVars',{'RT', 'soa', 'acc'});
    stats20.soa   = grpstats(d20,'soa',  {'mean', 'std'},'DataVars',{'RT', 'acc'});
    disp('Grouping Variable (unimodal - V2 only): trlid')
    disp(stats20.trlid)
    disp('Grouping Variable (unimodal - V2 only): soa')
    disp(stats20.soa)

    % A1
    stats01.trlid = grpstats(d01,'trlid',{'mean', 'std'},'DataVars',{'RT', 'soa', 'acc'});
    stats01.soa   = grpstats(d01,'soa',  {'mean', 'std'},'DataVars',{'RT', 'acc'});
    disp('Grouping Variable (unimodal - A1 only): trlid')
    disp(stats01.trlid)
    disp('Grouping Variable (unimodal - A1 only): soa')
    disp(stats01.soa)
    % A2
    stats02.trlid = grpstats(d02,'trlid',{'mean', 'std'},'DataVars',{'RT', 'soa', 'acc'});
    stats02.soa   = grpstats(d02,'soa',  {'mean', 'std'},'DataVars',{'RT', 'acc'});
    disp('Grouping Variable (unimodal - A2 only): trlid')
    disp(stats02.trlid)
    disp('Grouping Variable (unimodal - A2 only): soa')
    disp(stats02.soa)
    
    % PLOTS
    disp('Plotting results...')

    %% plot unimodal results
    if any(ismember(d.trlid, [2 3])) % V only
        % RT (boxplots)
        rtplot = figure;
        subplot(2,2,1)                          % all
        boxplot(d.RT,d.soa)
        title('Response times (V1 and V2)')
        xlabel('soa'); ylabel('RT');
        subplot(2,2,2)                          % V1 
        boxplot(d10.RT,d10.soa)
        title('Response times (V1)')
        xlabel('soa'); ylabel('RT');
        subplot(2,2,4)                          % V2
        boxplot(d20.RT,d20.soa)
        title('Response times (V2)')
        xlabel('soa'); ylabel('RT');
        set(gcf,'color','w');
        % set(find(0, 'Type', 'title'), 'FontSize', 18)

        % ACC 
        % 2 sd error bars - assuming normpdf, just to get an idea
        accplot = figure;
        subplot(2,2,1)                          % all
        scatter(stats.soa.soa, stats.soa.mean_acc)
        errorbar(stats.soa.soa, stats.soa.mean_acc, stats.soa.std_acc)
        title('Accuracy (V1 and V2)')
        xlabel('soa'); ylabel('acc')
        subplot(2,2,2)                          % V1
        scatter(stats10.soa.soa, stats10.soa.mean_acc)
        errorbar(stats10.soa.soa, stats10.soa.mean_acc, stats10.soa.std_acc)
        title('Accuracy (V1)')
        xlabel('soa'); ylabel('acc')
        subplot(2,2,4)                          % V2
        scatter(stats20.soa.soa, stats20.soa.mean_acc)
        errorbar(stats20.soa.soa, stats20.soa.mean_acc, stats20.soa.std_acc)
        title('Accuracy (V2)')
        xlabel('soa'); ylabel('acc')
        set(gcf,'color','w'); 
    end % v only
    
    if any(ismember(d.trlid, [4 7])) % A only
        % RT (boxplots)
        rtplot = figure;
        subplot(2,2,1)                          % all
        boxplot(d.RT,d.soa)
        title('Response times (A1 and A2)')
        xlabel('soa'); ylabel('RT');
        subplot(2,2,2)                          % A1 
        boxplot(d01.RT,d01.soa)
        title('Response times (A1)')
        xlabel('soa'); ylabel('RT');
        subplot(2,2,4)                          % A2
        boxplot(d02.RT,d02.soa)
        title('Response times (A2)')
        xlabel('soa'); ylabel('RT');
        set(gcf,'color','w');
        % set(find(0, 'Type', 'title'), 'FontSize', 18)

        % ACC 
        % 2 sd error bars - assuming normpdf, just to get an idea
        accplot = figure;
        subplot(2,2,1)                          % all
        scatter(stats.soa.soa, stats.soa.mean_acc)
        errorbar(stats.soa.soa, stats.soa.mean_acc, stats.soa.std_acc)
        title('Accuracy (A1 and A2)')
        xlabel('soa'); ylabel('acc')
        subplot(2,2,2)                          % A1
        scatter(stats01.soa.soa, stats01.soa.mean_acc)
        errorbar(stats01.soa.soa, stats01.soa.mean_acc, stats01.soa.std_acc)
        title('Accuracy (A1)')
        xlabel('soa'); ylabel('acc')
        subplot(2,2,4)                          % A2
        scatter(stats02.soa.soa, stats02.soa.mean_acc)
        errorbar(stats02.soa.soa, stats02.soa.mean_acc, stats02.soa.std_acc)
        title('Accuracy (A2)')
        xlabel('soa'); ylabel('acc')
        set(gcf,'color','w'); 
    end % v only
    
    try
        % Signal detection parameters
        fprintf('\n\nCalculating sensitivity and bias (UNISENSORY):\n')
        fprintf('\nSignal=V1, Noise=V2:\n');
        sdm23 = dfi_SDM( d, 2, 3 ); disp(sdm23);% V1 = signal, V2 = noise
        fprintf('\nSignal=V2, Noise=V1:\n');
        sdm32 = dfi_SDM( d, 3, 2 ); disp(sdm32);% V2 = signal, V1 = noise
        fprintf('\nSignal=A1, Noise=A2:\n');
        sdm47 = dfi_SDM( d, 4, 7 ); disp(sdm47);% A1 = signal, A2 = noise
        fprintf('\nSignal=A2, Noise=A1:\n');
        sdm74 = dfi_SDM( d, 4, 7 ); disp(sdm74);% A1 = signal, A2 = noise
    catch
        fprintf('\n\nProblem with calculating signal detection parameters. Check Palamedes installation.\n');
    end
end % unimodal condition


%% Staircase
if strcmp(s.stair.incl, 'yes')
    figure('color', 'w', 'position', [50 50 1700 800])
    for ii = 1:s.stair.nChains
        subplot(1,s.stair.nChains,ii)
        plot(s.ud{ii}.FF.x / (1/120), 'g.-', 'linewidth', 3, 'markersize', 30); hold on
        plot(s.ud{ii}.Fus.x/ (1/120), 'b.-', 'linewidth', 3, 'markersize', 30);
        plot(s.ud{ii}.Fis.x/ (1/120), 'r.-', 'linewidth', 3, 'markersize', 30);
        set(findall(gcf,'-property','LineWidth'), 'LineWidth', 1.5)
        title(sprintf('Chain group %i', ii), 'fontsize', 24);
        xlabel('Trial number', 'fontsize', 20); ylabel('Flips', 'fontsize', 20);
        set(findall(gcf, '-property', 'XTickLabel'), 'fontsize', 18); 
        grid on
    end
end % stair 



end % function

% eof





























