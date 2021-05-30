function [ varargout ] = dfi_eeg_analysis( s, eeg )
%varargout = dfi_eeg_analysis(s, eeg)
%   this function takes a structure with specifications as input (e.g. how
%   long should the trials be, which frequencies to analyze, etc...
%   as well as a data structure 'eeg' as given by fieldtrip after artifact
%   rejection (so the data is clean). See example for s details.
%
% OUTPUT:
%   varargout returns, depending on analyses, in this order:
%   'eeg'   (i.e. input data, but modified, e.g. filtered or padded),
%   'spect' (i.e. data after power spectral analysis) 
%   'tfr'   (i.e. data after time frequency analysis)
%   'phase' (i.e. complex fourier transform with average phase and
%                 inter-trial coherence (itc)
%   'syn'   (i.e. a structure with the following fields:
%           .plv       , phase locking value
%           .label     , channel labels in .plv
%           .pls       , proportion surrogate PLVs higher than orig PLV
%           .plsChA    , channels from
%           .plsChB    , channels to
%           .plsMat    , matrix of the form 
%                        timePoints x numPerm x channelsFrom x channelsTo
%
% EXAMPLE(s):
% - Note that you need to have a fieldtrip data structure 'eeg' in wm
%
% % Trials
% s.trls    = 'all';
% s.keeptrls= 'no';
% s.toilim  = [];
% 
% % Space and time
% s.chan    = 'all'; %{'O1'};                                            
% s.toi     = [-0.6 0.6];          % time window (Cecere: -0.5 0.5 s)
% s.padwin  = [-1.2 1.2];
% 
% % Frequency
% s.foi     = [1 35];              % frequency window to be plotted and bp window
% s.taper   = 'hanning';
% s.tapsmfr = 2;
% s.overlap = 0.01;                % TFR window overlap (smoothness) 
% s.t_ftwin = 1;
% s.pad     = 8;                   % zero padding (top up with 0 until pad s)
% s.filwin  = s.foi;               % bandpass window (Cecere: 3 40 Hz)
% s.filt    = 'no';
% s.filtord = 5;                   % filter order (always use two-pass to avoid phase shifts)
% s.filtplot= 'no';                % plot filter response (for each trial!)
% s.filtyp  = 'but';               % butterworth IIR filter
% s.alpha   = [8 14];              % alpha bounds 
% s.tfrbsltyp  = 'db';                % 'absolute', 'relative', 'relchange' (default = 'absolute')
% s.smoowin = 'uniform';           % or 'increase_with_frequency'
% s.trns_wdth = 0.15;              % Fslide filter transition width (default = 0.15 from Cohen)
% 
% 
% % Synchrony
% s.filtrng = [8 14];              % bandpass filter range (foi)
% s.filtcyc = 5;                   % cycles of foi the FIR should last
%                                  % The IR of an Nth-order discrete-time FIR filter 
%                                  % lasts exactly N + 1 samples (Wikipedia)
% % Analyses / plots
% s.psd     = 1;                   % - PSD (var name 'spect')
% s.psdsing = 1;                   %   in a single plot
% s.powscale= 'log';               %   'log', 'logsub', 'logperc' (default='powspctrm')
% s.psdmult = 0;                   %   or over scalp
% s.psdtopo = 0;                   %   or make topoplot
% 
% s.tfr     = 1;                   % - TFR (var name 'tfr')
% s.tfrsing = 1;                   %   in a single plot
% s.tfrmult = 0;                   %   or over scalp
% 
% s.phase   = 1;                   % - PHASE (var name 'phase')
% s.phasesp = 1;                   %   single plots
% s.phasemp = 0;                   %   arrange over scalp
% 
% s.plv     = 1;                   % - PLV (var name 'plv')
% s.pls     = 1;                   %   get phase locking statistics as well
% s.plsnperm= 200;                 %   number of permutations for significance test
% s.plvplot = 1;                   %   plot phase locking values for specified channels
% s.plsplot = 1;                   %   plot statistics
% s.plvchA  = {'O1', 'O2', 'Oz'};
% s.plvchB  = {'Cz'};
% 
% s.fslide  = 0;                   % - Frequency sliding (var name 'fslide')
% 
% [~, pow, tfr, phase, pl]  = dfi_eeg_analysis(s, eeg);
%
%
% // Steffen Buergers, sbuergers@gmail.com, July 2015
%

argcount = 1; % how many variables do I return?


%% 1.) --- TRIAL SELECTION ---
if s.redeftrl
    cfg        = []; 
    cfg.toilim = s.padwin; 
    cfg.trials = s.trls;
    if ~isempty(s.toilim) % redefine trial length (default = no)?
        cfg.toilim = s.toilim;
    end

    eeg = ft_redefinetrial(cfg, eeg);
end


%% 2.) --- FILTERING --- 
if ~strcmp(s.filt, 'no')
    cfg              = [];
    cfg.feedback     = s.feedback;
    cfg.bpfilter     = s.filt;
    cfg.bpfilttype   = s.filtyp;
    cfg.bpfreq       = s.filwin; 
    cfg.bpfiltord    = s.filtord;
    cfg.plotfiltresp = s.filtplot;

    eeg              = ft_preprocessing(cfg, eeg); 
end

varargout{argcount} = eeg;

if isempty(eeg.trial)
    fprintf('\n\nNo trials were selected, returning empty FT structure.....\n')
    return
end


%% 3.) --- FREQUENCY ANALYSIS: PSD ---
if s.psd
    
    cfg = [];
    cfg.feedback   = s.feedback;
    cfg.keeptrials = s.keeptrls;
    cfg.channel    = s.chan;
    cfg.output     = 'pow';
    cfg.method     = 'mtmfft';
    cfg.taper      = s.taper;
    cfg.tapsmofrq  = s.tapsmfr;
    cfg.foilim     = s.foi;
    cfg.pad        = s.pad;
    
    spect = ft_freqanalysis(cfg, eeg);
    
    % Powspctrm baseline corrections 
    % log:     logspectrum
    % logsub:  logspectrum minus average logspectrum
    % logperc: (logspectrum + average logspectrum) / average logspectrum * 100 (% deviation from mean)
    if strcmp(s.keeptrls, 'yes')
        [spect.log, spect.logsub, spect.logperc] = deal(zeros(size(spect.powspctrm)));
        for itrl = 1:size(spect.powspctrm,1)
            for ichan = 1:numel(spect.label)
                spect.log(itrl,ichan,:)     =  log10(spect.powspctrm(itrl,ichan,:));
                spect.logsub(itrl,ichan,:)  =  log10(spect.powspctrm(itrl,ichan,:)) - mean(log10(spect.powspctrm(itrl,ichan,min(spect.freq):max(spect.freq))),ndims(spect.powspctrm));
                spect.logperc(itrl,ichan,:) = (log10(spect.powspctrm(itrl,ichan,:)) + mean(log10(spect.powspctrm(itrl,ichan,min(spect.freq):max(spect.freq))),ndims(spect.powspctrm))) / ...
                    mean(log10(spect.powspctrm(itrl,ichan,min(spect.freq):max(spect.freq))),ndims(spect.powspctrm))*100;
            end
        end
    else
        for ichan = 1:numel(spect.label)
            spect.log(ichan,:)     =  log10(spect.powspctrm(ichan,:));
            spect.logsub(ichan,:)  =  log10(spect.powspctrm(ichan,:)) - mean(log10(spect.powspctrm(ichan,min(spect.freq):max(spect.freq))),2);
            spect.logperc(ichan,:) = (log10(spect.powspctrm(ichan,:)) + mean(log10(spect.powspctrm(ichan,min(spect.freq):max(spect.freq))),2))/mean(log10(spect.powspctrm(ichan,min(spect.freq):max(spect.freq))),2)*100;
        end
    end
    
    % PSD single plot
    if s.psdsing
        if ~isfield(s, 'powscale')
            s.powscale = 'powspctrm';
        end
        switch s.powscale
            case 'powspctrm', plbl = '�V^2';
            case 'log'      , plbl = 'log10(�V^2)';    
            case 'logsub'   , plbl = 'log10(�V^2)'; 
            case 'logperc'  , plbl = 'perc. change'; 
        end
        figure, hold on
        plot(spect.freq, spect.(s.powscale)(:, :)), legend(spect.label);
        plotspecs, shadearea(s.alpha);
        xlabel('Frequency (Hz)'), ylabel(plbl);
        hold off
    end

    % PSD scalp discribution / topoplot
    cfg.layout      = 'EEG1010.lay';
    cfg.parameter   = s.powscale;
    cfg.interactive = 'yes';
    cfg.colorbar    = 'yes';
    cfg.marker      = 'labels';
    if s.psdmult, figure, ft_multiplotER(cfg, spect); end; 
    if s.psdtopo, figure, ft_topoplotER( cfg, spect); end; 
    
    argcount = argcount + 1;
    varargout{argcount} = spect;
end % PSD


%% 4.) --- FREQUENCY ANALYSIS: Spectrogram ---
if s.tfr
    
    cfg              = [];
    cfg.feedback     = s.feedback;
    cfg.channel      = s.chan;
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = s.taper; % 1/fresol
    cfg.foi          = s.foi(1):(s.t_ftwin):s.foi(2);                         
    cfg.toi          = s.toi(1):s.overlap:s.toi(2);  
    cfg.pad          = s.pad;
    cfg.keeptrials   = s.keeptrls;
    if strcmp(s.taper, 'dpss')
        cfg.t_ftimwin = min(s.foi)./cfg.foi;   % (tw) time window decreases with frequency
        cfg.tapsmofrq = 5./cfg.foi;            % (fw) smoothing increases with frequency
    else
        if strcmp(s.smoowin, 'increase_with_frequency')
            cfg.t_ftimwin = s.tfr_scale_fact./cfg.foi; % time window decreases with frequency
                                                       % more smoothing for higher frequencies
        else
            cfg.t_ftimwin = repmat(s.t_ftwin, 1, numel(cfg.foi));           
        end
    end
    
    tfr = ft_freqanalysis(cfg, eeg);

    % TFR single plots
    if s.tfrsing
        cfg = [];
        cfg.maskstyle    = 'saturation';
        cfg.baseline     = s.tfrbsl;
        cfg.baselinetype = s.tfrbsltyp;
        cfg.colorbar     = 'no';
        figure; n = numel(tfr.label);
        for ich = 1:n
            if n < 3 && n > 1
                subplot(n,1,ich); 
            else 
                subplot(ceil(n/2),2,ich);
            end
            cfg.channel = tfr.label(ich);
            ft_singleplotTFR(cfg, tfr);
            title(tfr.label(ich)), xlabel('Time (s)'), ylabel('Frequency (Hz)'); 
            cbh   = colorbar; 
            if strcmp(s.tfrbsl, 'yes'), title(cbh, s.tfrbsltyp); end;
            plotspecs;         
            yt    = get(cbh, 'YTick');
            YTick = -max(abs(yt)):(yt(2) - yt(1)):max(abs(yt));
            caxis([min(YTick), max(YTick)]); set(cbh, 'YTick', YTick);
        end
    end
    
    % TFR scalp discribution
    if s.tfrmult
        cfg.layout       = 'EEG1010.lay';
        cfg.interactive  = 'yes';
        cfg.colorbar     = 'yes';
        cfg.marker       = 'labels';
        cfg.baseline     = s.tfrbsl;
        cfg.baselinetype = s.tfrbsltyp;
        figure, ft_multiplotTFR(cfg, tfr); 
    end; 
    
    argcount = argcount + 1;
    varargout{argcount} = tfr;
end % TFR


%% 5.) --- FREQUENCY ANALYSIS: Phase ---
if s.phase

    % get complex fourier spectrogram
    cfg            = [];
    cfg.feedback   = s.feedback;
    cfg.channel    = s.chan;
    cfg.output     = 'fourier';
    cfg.method     = 'mtmconvol'; 
    cfg.taper      = s.taper; % 1/fresol
    cfg.foi        = s.foi(1):(1/s.pad):s.foi(2);                         
    cfg.toi        = s.toi(1):s.overlap:s.toi(2);  
    cfg.pad        = s.pad;
    cfg.keeptrials = 'yes';
    if strcmp(s.taper, 'dpss')
        cfg.t_ftimwin = min(s.foi)./cfg.foi;   % (tw) time window decreases with frequency
        cfg.tapsmofrq = 1/5 * cfg.foi;         % (fw) smoothing increases with frequency
    else
        if strcmp(s.smoowin, 'increase_with_frequency')
            cfg.t_ftimwin = s.tfr_scale_fact./cfg.foi; % time window decreases with frequency
                                                       % more smoothing for higher frequencies
        else
            cfg.t_ftimwin = repmat(s.t_ftwin, 1, numel(cfg.foi));           
        end
    end
    
    stft = ft_freqanalysis(cfg, eeg);
    
    % get itc and average phase
    n = numel(stft.label);
    if s.phasesp, figure; end;
    [stft.avgph, stft.itc] = getITC(stft);
    
    % PHASE single plots
    if s.phasesp
        for ich = 1:n
            if n < 3 && n > 1
                subplot(n,1,ich);
            else
                subplot(ceil(n/2),2,ich);
            end
            clim = s.toi;
            if n == 1
                imagesc( stft.time, stft.freq, stft.itc, clim );
            else
                imagesc( stft.time, stft.freq, squeeze(stft.itc(1,:,:)), clim );
            end
            title(stft.label(ich)), xlabel('Time (s)'), ylabel('Frequency (Hz)');
            colorbar; plotspecs
        end
    end
    
%     % Use FieldTrip example code to get ITC as well (includes linear ITC)
%     % make a new FieldTrip-style data structure containing the ITC
%     % copy the descriptive fields over from the frequency decomposition
%     itc = [];
%     itc.label     = stft.label;
%     itc.freq      = stft.freq;
%     itc.time      = stft.time;
%     itc.dimord    = 'chan_freq_time';
% 
%     F = stft.fourierspctrm;   % copy the Fourier spectrum
%     N = size(F,1);            % number of trials
% 
%     % compute inter-trial phase coherence (itpc) 
%     itc.itpc      = sum(F,1) ./ sum(abs(F),1);
%     itc.itpc      = abs(itc.itpc);     % take the absolute value, i.e. ignore phase
%     itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
% 
%     % compute inter-trial linear coherence (itlc)
%     itc.itlc      = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
%     itc.itlc      = abs(itc.itlc);     % take the absolute value, i.e. ignore phase
%     itc.itlc      = squeeze(itc.itlc); % remove the first singleton dimension
%     
%     stft.itcft.itpc = itc.itpc;
%     stft.itcft.itlc = itc.itlc;
    
    % PHASE scalp discribution
    if s.phasemp
        figure;
        dumm = stft;
        dumm.dimord = 'chan_freq_time';

        % plot at channel location
        cfg = [];
        cfg.zparam   = 'itc';
        cfg.colorbar = 'yes';
        cfg.layout   = 'EEG1010.lay';
        ft_multiplotTFR(cfg, dumm)
    end
    
    argcount = argcount + 1;
    varargout{argcount} = stft;
end % PHASE    

    
%% 6.) --- FREQUENCY ANALYSIS: Synchrony ---
if s.plv
% NOTE:
% As you have probably noticed in the plot from the above example, the PLV 
% between two random signals is spuriously large in the first 100 ms. While 
% using FIR filtering and/or hilbert transform, it is good practice to 
% discard both ends of the signal (same number of samples as the order of 
% the FIR filter, or more).
% 
% Also note that in order to extract the PLV between channels 17 and 20, 
% use plv(:, 17, 20, :) and NOT plv(:, 20, 17, :). The smaller channel 
% number is to be used first.

    % select relevant channels
    if exist('s.plvchA', 'var') && exist('s.plvchB', 'var')
        cfg            = [];
        cfg.channel    = [s.plvchA, s.plvchB];

        eegplv         = ft_preprocessing(cfg, eeg); 
    else
        eegplv         = eeg;
    end

    % The IR of an Nth-order discrete-time FIR filter lasts exactly N+1 samples (Wikipedia)
    Fs         = eeg.fsample;
    avgfoi     = mean(s.filtrng);
    filtord    = round(1/avgfoi * Fs * s.filtcyc);
    
    % select data segment of interest 
    cfg        = [];
    cfg.toilim = [s.toi(1) - filtord/Fs, s.toi(2) + filtord/Fs]; 
    if cfg.toilim(1) < min(eeg.time{1}) || cfg.toilim(2) > max(eeg.time{1})
        error('PLV: Your data segment is not large enough for adequate PLV estimation!')
    end
    eegplv     = ft_redefinetrial(cfg, eegplv);
    
    eegData    = zeros([size(eegplv.trial{1}), numel(eegplv.trial)]);
    for i = 1:numel(eegplv.trial)
        try
            eegData(:,:,i) = eegplv.trial{i}(:,1:end);
        catch
            eegData(:,:,i) = eegplv.trial{i}(:,1:end-1);
        end
    end
    srate          = eegplv.fsample; 
    filtSpec.order = filtord;
    filtSpec.range = s.filtrng; 

    % calculate PLV                        
    % the output is of the form     'plv'    'chanFrom'    'chanTo'
    [plval] = pn_eegPLV(eegData, srate, filtSpec);
    
    % Decide what constitutes a high PLV by computing surrogate data
    if s.pls
        if numel(s.plvchB) * numel(s.plvchA) > 20
            fprintf('\n\nWARNING: Omitting calculation of PLV statistics to avoid memory saturation. \nProvide fewer channels pairs per call (max = 20)...\n\n');
            surplvMat = [];
            s.pls = 0;
        else
            if ~isfield(s, 'plsnperm'), error('please specify s.nperm for plotting PLS.'), end;
            surplvMat = zeros([size(squeeze(eegData(1,:,1)),2), s.plsnperm, numel(s.plvchA), numel(s.plvchB)]);
            for iFrom = 1:numel(s.plvchA)
                for iTo = 1:numel(s.plvchB)
                    elA = strcmp(s.plvchA(iFrom), eegplv.label);
                    elB = strcmp(s.plvchB(iTo), eegplv.label);
                    surplvMat(:,:,iFrom,iTo) = getPLS( squeeze(eegData(elA,:,:)), squeeze(eegData(elB,:,:)), ...
                                                       srate, filtSpec, s.plsnperm );
                end
            end
        end
    end
    
    if s.pls
        % the proportion of surrogate samples with higher PLVs constitutes our
        % measure of significance
        pls = zeros([size(plval,1), numel(s.plvchA), numel(s.plvchB)]);
        for iA = 1:numel(s.plvchA)
            for iB = 1:numel(s.plvchB)
                idxA = find(strcmp(eegplv.label, s.plvchA(iA)));
                idxB = find(strcmp(eegplv.label, s.plvchB(iB)));
                if idxA < idxB
                    numHigher = sum( surplvMat(:,:,iA,iB) >= repmat(plval(:,idxA,idxB),[1 s.plsnperm]) ,2);
                else
                    numHigher = sum( surplvMat(:,:,iA,iB) >= repmat(plval(:,idxB,idxA),[1 s.plsnperm]) ,2);
                end
                pls = numHigher / s.plsnperm;
            end
        end
    end

    % PLOT: PLV
    if s.plvplot
        chA = find(ismember(eegplv.label, s.plvchA));
        chB = find(ismember(eegplv.label, s.plvchB));
        figure; i = 0;
        for iA = 1:numel(chA)
            for iB = 1:numel(chB)
                i = i + 1;
                subplot(numel(chA), numel(chB), i)
                if chA(iA) < chB(iB)
                    plvData = squeeze(plval(:, chA(iA), chB(iB)));
                else
                    plvData = squeeze(plval(:, chB(iB), chA(iA)));
                end
                plot((s.toi(1)*Fs:s.toi(2)*Fs-1)/srate, ...
                      plvData(floor(numel(plvData)/2+s.toi(1)*Fs) : floor((numel(plvData)-1)/2+s.toi(2)*Fs-1)));
                xlabel('Time (s)'); ylabel('PLV'); xlim(s.toi);
                title(sprintf('%s - %s', eegplv.label{chA(iA)}, eegplv.label{chB(iB)}));
            end
        end; plotspecs
    end
    
    % PLOT: PLV and surrogate distribution
    if s.plsplot
        chA = find(ismember(eegplv.label, s.plvchA));
        chB = find(ismember(eegplv.label, s.plvchB));
        figure; i = 0;
        for iA = 1:numel(chA)
            for iB = 1:numel(chB)
                i = i + 1;
                subplot(numel(chA), numel(chB), i)
                if chA(iA) < chB(iB)
                    plvData = squeeze(plval(:, chA(iA), chB(iB)));
                else
                    plvData = squeeze(plval(:, chB(iB), chA(iA)));
                end
                hold on
                plot((s.toi(1)*Fs:s.toi(2)*Fs-1)/srate, ...
                      surplvMat(floor(numel(plvData)/2+s.toi(1)*Fs) : floor((numel(plvData)-1)/2+s.toi(2)*Fs-1), :, iA, iB), 'color', [.5 .5 .5]);
                plot((s.toi(1)*Fs:s.toi(2)*Fs-1)/srate, ...
                      plvData(  floor(numel(plvData)/2+s.toi(1)*Fs) : floor((numel(plvData)-1)/2+s.toi(2)*Fs-1)));
                xlabel('Time (s)'); ylabel('PLV');  xlim(s.toi);
                title(sprintf('%s - %s', eegplv.label{chA(iA)}, eegplv.label{chB(iB)}));
                hold off
            end
        end; plotspecs
    end
    
    % collect output
    syn.plv    = plval;
    syn.label  = eegplv.label;
    if exist('pls', 'var'),       syn.pls = pls;          end;
    if exist('s.plvchA', 'var'),  syn.plsChA = s.plvchA;  end;
    if exist('s.plvchB', 'var'),  syn.plsChB = s.plvchB;  end;
    if exist('surplvMat', 'var'), syn.plsMat = surplvMat; end;
    
    argcount = argcount + 1;
    varargout{argcount} = syn;
end % SYNCHRONY


%% 7.) --- FREQUENCY ANALYSIS: Frequency sliding ---
if s.fslide
% Following the method outlined in Cohen (2014):
%     Fluctuations in oscillation frequency control spike timing and coordinate 
%     neural networks. The Journal of Neuroscience, 2 July 2014, 34(27): 8988-8998; 
%     doi: 10.1523/JNEUROSCI.0261-14.2014
%
% Frequency sliding occurs from single neurons (f-I-curve) to large
% ensembles and slower oscillations lead to a lower spiking threshold and
% more spike timing variability. In large scale EEG it can index functional
% connectivity between regions. 
%
% Instantaneous frequency is defined as the change in phase per unit time, 
% i.e. the temporal derivative of the phase angle time-series. 

    % Use fast median filter (from file exchange)
    % addpath('E:\nth_element_0.87')
    fprintf('\nFrequency sliding: Filtering data with tukeywin=[%i %i].....\n', ...
            s.alpha(1), s.alpha(2)); tic
    
    if ~strcmp(s.chan, 'all')
        cfg            = [];
        cfg.channel    = s.chan;
        eeg            = ft_preprocessing(cfg, eeg);
    end

    % Step 1:     Get raw EEG data (already happened)
    % Step 2 & 3: Plateau shaped band pass cosine filter (I use a low and high
    % pass filter as no pleateau filter is available in fieldtrip)
    if strcmp(s.tukeyfilt, 'firws')
        cfg               = [];
        cfg.feedback      = s.feedback;
        cfg.lpfilter      = 'yes';
        cfg.lpfilttype    = 'firws';
        cfg.lpfiltwintype = 'hann';
        cfg.lpfiltord     = 512;%768;%512;%256;%
        cfg.lpfreq        = s.alpha(2); 
        cfg.hpfilter      = 'yes';
        cfg.hpfilttype    = 'firws';
        cfg.hpfiltwintype = 'hann';
        cfg.hpfiltord     = 512;%768;%512;%256;%
        cfg.hpfreq        = s.alpha(1);
        cfg.plotfiltresp  = 'no';

        eeg2              = ft_preprocessing(cfg, eeg);

        % Plot filter properties
        if s.fsplot 
            cfg.trials = 1;
            cfg.plotfiltresp = 'yes';
            ft_preprocessing(cfg, eeg2);
        end
        eeg3 = eeg2;
    
    elseif strcmp(s.tukeyfilt, 'firls')
        % filter data (code MXCohen 2014)
        % create some useful variables that are otherwise in the eeglab struct
        fs    = eeg.fsample;
        nch   = numel(eeg.label); 
        ntrls = numel(eeg.trial);
        ntpts = numel(eeg.time{1});

        % transform my fieldtrip data to eeglab format (i.e. a matrix of ch x tp x trl)
        eegdata = zeros([nch, ntpts, ntrls]);
        for itrl = 1:numel(eeg.trial)
            eegdata(:,:,itrl) = eeg.trial{itrl};
        end
        % apply a band-pass filter with 15% transition zones.
        trns_wdth   = s.trns_wdth; 
        idealrsp    = [ 0 0 1 1 0 0 ];
        filtfrqbnds = [ 0 (1-trns_wdth)*s.alpha(1) s.alpha(1) s.alpha(2) s.alpha(2)*(1+trns_wdth) fs/2 ]/(fs/2);
        filt_ord    = round(3*(fs/s.alpha(1)));
        filtwghts   = firls(filt_ord,filtfrqbnds,idealrsp);
        
        
%         % plot filter response
%         nyquist = fs/2;
%         figure('color', 'w')
%         plot(filtwghts,'r')
%         xlabel('Time')
%         title('Filter response'); plotspecs
%         
%         figure('color', 'w')
%         plot(filtfrqbnds*nyquist,idealrsp,'r'); hold on
%         fft_filtkern  = abs(fft(filtwghts));
%         fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
%         hz_filtkern   = linspace(0,nyquist,ceil(length(fft_filtkern)/2));
%         plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'b')
%         set(gca,'ylim',[-.1 1.1],'xlim',[0 nyquist])
%         xlabel('Frequency');
%         legend({'ideal';'best fit'}); plotspecs; xlim([1 35])
%         

        % this part does the actual filtering
        filterdata = zeros(size(eegdata));
        for chani = 1:nch
            filterdata(chani,:,:) = reshape( filtfilt(filtwghts,1,double(reshape(eegdata(chani,:,:),1,ntpts*ntrls))) ,ntpts, ntrls);
        end
        
        % adjust to other method's flow
        eeg3 = eeg;

        % reformat to fieldtrip structure
        for itrl = 1:numel(eeg.trial)
            eeg3.trial{itrl} = filterdata(:,:,itrl);
        end
    end % filter selection
    
    fprintf('Frequency sliding: Filtering took %d seconds.....\n', round(toc)); 
    fprintf('Frequency sliding: Applying hilbert transform.....\n'); 
    % Step 4 & 5: apply hilbert transform and get phase angle time series
    option = 'complex';
    for itrl = 1:numel(eeg3.trial)
        % ft does exactly this: dat = transpose(hilbert(transpose(dat)));
        eeg3.hilb{itrl}  = ft_preproc_hilbert(eeg3.trial{itrl}, option);
        for ich = 1:length(eeg3.label)
            eeg3.phase{itrl}(ich,:) = angle(eeg3.hilb{itrl}(ich,:));
        end
    end
    
    fprintf('Frequency sliding: Get temporal derivative of phase angle series.....\n'); 
    % Step 6: % Get temporal derivative of phase-angle time series and apply
    % median filter to quench artifacts
    % Convert temporal derivative of phase angle time-series to Hertz
    % Hz(t) = s*(Phi(t) - Phi(t-1))/2pi   
    for itrl = 1:numel(eeg3.trial)
        % initialize new structure for freq sliding data
        dimen = size(eeg3.trial{itrl});
        eeg3.fslide_unfilt{itrl} = NaN([dimen(1), dimen(2)-1]);
        % create new time axis (we loose a point due to the diff function)
        eeg3.t_unfilt{itrl} = eeg3.time{itrl}(1:end-1);
        % get instantaneous frequency
        eeg3.fslide_unfilt{itrl} = diff(eeg3.fsample*unwrap(eeg3.phase{itrl},[],1),1,2)/(2*pi);
    end
    
    % Implementation by Cohen, 2014
    % Apply 10 median filters with different window lengths and then take 
    % the median of the median filters and average over trials
    n_order = 10;
    filtord = round(linspace(10,400,n_order))/2;
    %filtord = round(linspace(10,400,n_order)); % TRY MEDFILT OF JS
    filtord = round( filtord/(1000/eeg.fsample) );
    timeoi  = (s.padwin(1)*1000:(1000/eeg.fsample):s.padwin(2)*1000)/1000;
    timeoiidx = dsearchn(eeg3.time{1}',timeoi');
    freqslide = zeros(length(eeg3.label),length(timeoi));
    %fslidemed = zeros(numel(eeg3.label),length(filtord),length(timeoi),numel(eeg3.trial)); % TRY MEDFILT OF JS
    frslide   = zeros(numel(eeg3.label), numel(eeg3.fslide_unfilt{1}(1,:)), numel(eeg3.trial));
    for itrl = 1:numel(eeg3.trial)
        frslide(:,:,itrl) = eeg3.fslide_unfilt{itrl}; %fslide size:  64        1229         330
    end
    
    
    
    
    % manual median filter from Cohen (almost identical to JS)
    fprintf('Frequency sliding: Applying median filters.....\n'); tic
    % median filter (only on requested time points to decrease computation time)
    for oi=1:n_order
        for ti=1:length(timeoi)
%             % use compiled fast_median if available
%             for triali=1:numel(eeg3.trial)
%                 fslidemed(:,oi,ti,triali) = fast_median(frslide(:,max(timeoiidx(ti)-filtord(oi),1):min(timeoiidx(ti)+filtord(oi),numel(eeg.time{1})-1),triali)');
%             end
            % use 'manual median' otherwise
            temp = sort(frslide(:,max(timeoiidx(ti)-filtord(oi),1):min(timeoiidx(ti)+filtord(oi),numel(eeg.time{1})-1),:),2);
            fslidemed(:,oi,ti,:) = temp(:,floor(size(temp,2)/2)+1,:);
        end
    end
    
    % Show that it is fine to have a larger padding window than actual
    % data for my small window of interest (-0.6 to -0.1):
%     figure
%     t_eegwin = eeg3.time{1}';
%     plot(timeoi, ones(length(timeoi),1)); hold on
%     plot(t_eegwin, ones(length(eeg3.time{1}'),1), 'r');
%     plot(t_eegwin(t_eegwin >= -0.6 & t_eegwin <= -0.1), ones(sum(t_eegwin >= -0.6 & t_eegwin <= -0.1),1), 'c')
%     plot(timeoi, squeeze(fslidemed(1,1,:,1)), 'k')
    
%     % median filter, exactly similar to Samaha's script % TRY MEDFILT OF JS
%     fslidemed = zeros(numel(eeg3.label), length(filtord),size(frslide,2),numel(eeg3.trial));
%     for oi=1:n_order
% 
%         fslidemed(:,oi,1:size(frslide,2),:) = medfilt1(frslide,filtord(oi),[],2);
%         
%     end

    % the final step is to take the mean of medians
    fprintf('Frequency sliding: Get median of the %i median filters and average over trials.....\n', n_order);
    medmatrix      = median(fslidemed,2);
    freqslide      = nanmean(medmatrix,4);
    
    fprintf('Frequency sliding: Median filtering took %d seconds.....\n\n\n', round(toc));
    % collect output
    eeg3.fslide_med    = squeeze(medmatrix);
    eeg3.fslide_avg    = freqslide;
    eeg3.fslide_toi    = timeoi';
    eeg3.fslide_toiidx = timeoiidx;
    
    % Plot processing steps and final outcome for example trial and channel
    if s.fsplot 
        figure('color','w');
        t = eeg3.time{1};
        subplot(611); plot(eeg.time{1}, eeg.trial{1}(1,:)); title('Step 1: Raw signal'); set(gca,'xlim',s.toi); set(gca, 'XtickLabel', []); xlim([-0.5 0]); grid on
        subplot(612); plot(t, eeg3.trial{1}(1,:)); title('Step 2: Filtered data, (8-14 Hz)'); set(gca,'xlim',s.toi); set(gca, 'XtickLabel', []); xlim([-0.5 0]); grid on
        subplot(613); plot(t, eeg3.phase{1}(1,:)); title('Step 3: Phase'); set(gca,'xlim',s.toi); set(gca, 'XtickLabel', []); xlim([-0.5 0]); grid on
        subplot(614); plot(t(1:end-1), eeg3.fslide_unfilt{1}(1,:)); title('Step 4: Instantaneous frequency'); set(gca,'xlim',s.toi, 'ylim', [s.alpha(1)-3, s.alpha(2)+3]); set(gca, 'XtickLabel', []); xlim([-0.5 0]); grid on
%         subplot(615); plot(t(1:end-1), squeeze(medmatrix(1,:,:,1))); title('Step 5: Median filtered'); set(gca,'xlim',s.toi); set(gca, 'XtickLabel', []); xlim([-0.5 0]); grid on
%         subplot(616); plot(t(1:end-1), freqslide(1,:)); xlabel('time'); title('Step 6: Averaged over trials'); set(gca,'xlim',s.toi); xlim([-0.5 0]); grid on
%        
        subplot(615); plot(eeg3.fslide_toi, squeeze(medmatrix(1,:,:,1))); title('Step 5: Median filtered'); set(gca,'xlim',s.toi); set(gca, 'XtickLabel', []); xlim([-0.5 0]); grid on
        subplot(616); plot(eeg3.fslide_toi, freqslide(1,:)); xlabel('time'); title('Step 6: Averaged over trials'); set(gca,'xlim',s.toi); xlim([-0.5 0]); grid on
        suptitle(sprintf('Trial %d, Channel %s', 1, eeg3.label{1}));
    end
    
    argcount = argcount + 1;
    varargout{argcount} = eeg3;
    
end % Frequency sliding


return

% eof


    




