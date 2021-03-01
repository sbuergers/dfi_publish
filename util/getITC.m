function [ avg_phase, itc ] = getITC( ft_data, toi, foi, coi, trloi )
%[ avg_phase, itc ] = getITC( ft_data, toi, foi, coi, trloi )
% takes a data structure from the output of ft_freqanalysis with 
% output='fourier' and calculates the inter-trial coherence and average
% phase.
%
% As the input can have up to four dimensions (trl, chan, freq, time), the 
% output will correspondingly be a 3d matrix, with the trial dimension
% removed. Optionally 'toi', 'foi' and 'coi' can be used to index the time
% point when phase coherence should be estimated, as well as frequencies
% and channels of interest to be returned.
%
% Note that in a complex vector in stft: -1.6734  +  0.0837i
%                                          ampl.    phase(rad)
%                                         |       .
%                                         |      .            
%                                         |     .
%                                         |   C.
%                                        B|   . ampl. = C, pow = C^2 (abs)
%                                         |  .
%                                         | .
%                                         |. ? = phase angle (cmplx)
%                                         .________
%                                             A
%
% ITC explanation:
% intertrial coherence (ITC, also called phase-locking factor), is a measure 
% that tells you how often a certain phase occurs at a certain time-period. 
% If at time=x, the phase always is y, then the ITC will be 1. If y is always 
% different at time x, then the ITC at that time-point is 0. In a complex plane 
% all datapoints are vectors with a certain phase (the angle) and a certain 
% amplitude. If we normalize all amplitudes (set them to 1), then all 
% datapoints will be on the unit-circle. If we then compute the (complex) 
% mean, we will find a vector with amplitude somewhere in between 0 and 1, 
% with a phase that represents the average phase at that time-point. The 
% amplitude of this vector is the ITC. 
%
% ITC explanation and code adapted from David Meijer, 2015
%

% check input data
if ~isfield(ft_data, 'fourierspctrm')
    error('Wrong input data format. Please provide field ´fourierspctrm´');
end;

% augment readability (stft = short time fourier transform)
stft = ft_data.fourierspctrm;

% select time point, frequencies and channels of interest
if exist('trloi', 'var'), 
    stft = stft(trloi,:,:,:); 
    numtrls = numel(trloi);
else
    numtrls = size(ft_data.trialinfo,1);
end;
if exist('toi', 'var'), stft = stft(:, :, :, toi); end;
if exist('foi', 'var'), stft = stft(:, :, foi, :); end;
if exist('coi', 'var'), 
    stft = stft(:, coi, :, :); 
    numch = numel(coi); 
else
    numch = numel(ft_data.label);
end;

% go through the stft and normalize the amplitudes (C = 1)
disp('Calculating inter-trial coherence...')
for c = 1:numch
    disp(['Channel: ', ft_data.label{c}])
    for i = 1:numtrls
          stft_norm(i,c,:,:) = stft(i,c,:,:)./abs(stft(i,c,:,:));
    end
end

% get modulus and phase angle from average of normalized amplitude vectors
itc       = squeeze(   abs(mean(stft_norm,1))  );
avg_phase = squeeze( angle(mean(stft_norm,1))  );

return 
% eof



















