function [s] = dfi_eeg_defaults
% [s] = dfi_eeg_defaults returns the structure s including defaults for the
% function dfi_eeg_analysis. It is recommended to always call this function
% first, make changes to the defaults in s and then call dfi_eeg_analysis.

% redefine trial at the beginning?
s.redeftrl = 1;

% Trials
s.trls    = 'all';
s.keeptrls= 'no';
s.toilim  = [];

% Verbosity of FieldTrip
s.feedback= 'yes';

% Space and time
s.chan    = 'all'; %{'O1'};                                            
s.toi     = [-0.6 0.6];          % time window (Cecere: -0.5 0.5 s)
s.padwin  = [-1.2 1.2];

% Frequency
s.foi     = [1 35];              % frequency window to be plotted and bp window
s.taper   = 'hanning';
s.tapsmfr = 2;
s.overlap = 0.01;                % TFR window overlap (smoothness) 
s.t_ftwin = 1;
s.pad     = 8;                   % zero padding (top up with 0 until pad s)
s.filwin  = s.foi;               % bandpass window (Cecere: 3 40 Hz)
s.filt    = 'no';
s.filtord = 5;                   % filter order (always use two-pass to avoid phase shifts)
s.filtplot= 'no';                % plot filter response (for each trial!)
s.filtyp  = 'but';               % butterworth IIR filter
s.alpha   = [8 14];              % alpha bounds 
s.tfrbsl  = 'yes';
s.tfrbsltyp  = 'db';                % 'absolute', 'relative', 'relchange' (default = 'absolute')
s.smoowin = 'uniform';           % or 'increase_with_frequency'
s.tfr_scale_fact = 5;            % for cfg.t_ftimwin = s.tfr_scale_fact./cfg.foi;
s.trns_wdth = 0.15;              % Fslide filter transition width (default = 0.15 from Cohen)



% Synchrony
s.filtrng = [8 14];              % bandpass filter range (foi)
s.filtcyc = 5;                   % cycles of foi the FIR should last
                                 % The IR of an Nth-order discrete-time FIR filter 
                                 % lasts exactly N + 1 samples (Wikipedia)
% Analyses / plots
s.psd     = 1;                   % - PSD (var name 'spect')
s.psdsing = 0;                   %   in a single plot
s.powscale= 'log';               %   'log', 'logsub', 'logperc' (default='powspctrm')
s.psdmult = 0;                   %   or over scalp
s.psdtopo = 0;                   %   or make topoplot

s.tfr     = 1;                   % - TFR (var name 'tfr')
s.tfrsing = 0;                   %   in a single plot
s.tfrmult = 0;                   %   or over scalp

s.phase   = 1;                   % - PHASE (var name 'phase')
s.phasesp = 0;                   %   single plots
s.phasemp = 0;                   %   arrange over scalp

s.plv     = 1;                   % - PLV (var name 'plv')
s.pls     = 1;                   %   get phase locking statistics as well
s.plsnperm= 200;                 %   number of permutations for significance test
s.plvplot = 1;                   %   plot phase locking values for specified channels
s.plsplot = 1;                   %   plot statistics
s.plvchA  = {'O1', 'O2', 'Oz'};
s.plvchB  = {'Cz'};

s.fslide  = 0;                   % - Frequency sliding (var name 'fslide')

return

% eof







