function [ surPLVmat ] = getPLS( dataChA, dataChB, srate, filtSpec, nperm)
%surPLVmat = getPLS( dataChA, dataChB, srate, filtSpec, selectArr );
% create surrogate phase locking values and compare to original phase
% locking values for statistical evaluation for channels A and B.
%
% As the distribution of phase locking values calculated from the data is
% not known (we cannot assume that it is uniform for example), we are going
% to create a surrogate distribution that we can use for statistical
% comparison. 
%
% INPUT:
%   dataChA / dataChB are 2D matrices with numTimePoints x nTrl
%   srate is the sampling rate of the EEG data
%   filtSpec is the filter specification to filter the EEG signal in the
%     desired frequency band of interest. It is a structure with two
%     fields, order and range. 
%       Range specifies the limits of the frequency
%     band, for example, put filtSpec.range = [35 45] for gamma band.
%       Specify the order of the FIR filter in filtSpec.order. A useful
%     rule of thumb can be to include about 4 to 5 cycles of the desired
%     signal. For example, filtSpec.order = 50 for eeg data sampled at
%     500 Hz corresponds to 100 ms and contains ~4 cycles of gamma band
%     (40 Hz).
%   nperm (default = 200) is the number of surrogate phase locking values
%
% OUTPUT:
%   surPLVmat is a numTimePoints x nPerm matrix containing surrogate phase
%     locking values
%
% H0: Signal A is independent from Signal B
% Intuitively the PLV is computed between two electrodes like such:
% 1.) at each time-point in signal A and signal B the instantaneous phase
%     is calculated and the difference between them is computed.
% 2.) This is done for all trials and finally an average difference phase
%     is computed. 
% 3.) to test whether this is significant, the trials in electrode B are
%     jumbled up randomly and steps 1. and 2. are repeated with that data.
% 4.) This is repeated 200 (or more) times, after which the number of repetitions on
%     which the surrogate PLV exceeded the original PLV is taken as an index of
%     how probable this PLV was even when there is no relationship between the
%     two signals (note that this assumes that there is variation between
%     trials and that the PLV does not stay constant throughout)
%
% Code adapted from:
% Praneeth Namburi
% http://uk.mathworks.com/matlabcentral/fileexchange/31600-phase-locking-value/content/pn_eegPLV.m
% Cognitive Neuroscience Lab, DUKE-NUS
% 01 Dec 2009
% 
% Present address: Neuroscience Graduate Program, MIT
% email:           praneeth@mit.edu
%-----------------------------------------
% Steffen Buergers, sbuergers@gmail.com
%

disp(['Calculating PLS with ' mat2str(nperm) ' iterations...']);

if ~exist('nperm', 'var')
    nperm = 200;
end

nPerm = nperm;
nTrl  = size(dataChA, 3);

filtPts  = fir1(filtSpec.order, 2/srate*filtSpec.range);
eegfiltA = filter(filtPts, 1, dataChA, [], 1);
eegfiltB = filter(filtPts, 1, dataChB, [], 1);

phaseA    = angle(hilbert(eegfiltA));
eegPhaseB = angle(hilbert(eegfiltB));

surplv = zeros(size(phaseA, 1), nPerm);
for iPerm = 1:nPerm
    phaseB        = eegPhaseB(:,randperm(size(eegPhaseB, 2)));
    surplv(:, iPerm) = abs(sum(exp(1i*(phaseA - phaseB)), 2)) / size(phaseA, 2);
end

surPLVmat = surplv;

return;


% eof


















