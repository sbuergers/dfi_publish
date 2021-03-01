function [plv] = pn_eegPLV(eegData, srate, filtSpec, selectArr)
% Computes the Phase Locking Value (PLV) for an EEG dataset.
%
% Input parameters:
%   eeg is a 3D matrix nCh x numTimePoints x nTrl
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
%   selectArr (OPTIONAL) is a logical 2D matrix of size - nTrl x
%     nCond. For example, if you have a 250 trials in your EEG
%     dataset and the first 125 correspond to the 'attend' condition and
%     the last 125 correspond to the 'ignore' condition, then use
%     selectArr = [[true(125, 1); false(125, 1)],...
%       [false(125, 1); true(125, 1)]];
%
% Output parameters:
%   plv is a 4D matrix - 
%     numTimePoints x nCh x nCh x nCond
%   If 'selectArr' is not specified, then it is assumed that there is
%   only one condition and all trials belong to that condition.
%
%--------------------------------------------------------------------------
% Example: Consider a 28 channel EEG data sampled @ 500 Hz with 231 trials,
% where each trial lasts for 2 seconds. You are required to plot the phase
% locking value in the gamma band between channels Fz (17) and Oz (20) for
% two conditions (say, attend and ignore). Below is an example of how to
% use this function.
%
%   eeg = rand(28, 1000, 231); 
%   srate = 500; %Hz
%   filtSpec.order = 50;
%   filtSpec.range = [35 45]; %Hz
%   selectArr = rand(231, 1) >= 0.5; % attend trials
%   selectArr(:, 2) = ~selectArr(:, 1); % ignore trials
%   [plv] = pn_eegPLV(eeg, srate, filtSpec, selectArr);
%   figure; plot((0:size(eeg, 2)-1)/srate, squeeze(plv(:, 17, 20, :)));
%   xlabel('Time (s)'); ylabel('Plase Locking Value');
%
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
%--------------------------------------------------------------------------
% 
% Reference:
%   Lachaux, J P, E Rodriguez, J Martinerie, and F J Varela. 
%   “Measuring phase synchrony in brain signals.” 
%   Human brain mapping 8, no. 4 (January 1999): 194-208. 
%   http://www.ncbi.nlm.nih.gov/pubmed/10619414.
% 
%--------------------------------------------------------------------------
% Written by: 
% Praneeth Namburi
% Cognitive Neuroscience Lab, DUKE-NUS
% 01 Dec 2009
% 
% Present address: Neuroscience Graduate Program, MIT
% email:           praneeth@mit.edu

nCh  = size(eegData, 1);
nTrl = size(eegData, 3);
if ~exist('selectArr', 'var')
    selectArr = true(nTrl, 1);
else
    if ~islogical(selectArr)
        error('Data selection array must be a logical');
    end
end
nCond = size(selectArr, 2);

disp('Filtering data...');
filtPts  = fir1(filtSpec.order, 2/srate*filtSpec.range);
% Filter frequency response
% freqz(filtPts,1,srate)
% filtPtsTest  = fir1(filtSpec.order, 2/1000*filtSpec.range);
% fvtool(filtPtsTest)
eegfilt  = filter(filtPts, 1, eegData, [], 2);
% compare to filtfilt function
% for ich = 1:size(eegData, 1)
%     for itrl = 1:size(eegData, 3)
%         eegfilt2(ich,:,itrl)  = filtfilt(filtPts, 1, eegData(ich,:,itrl));
%     end
% end
% figure;
% subplot(2,1,1); plot(eegfilt(ich, :, itrl));
% subplot(2,1,2); plot(eegfilt(ich, :, itrl));

eegPhase = eegfilt;

disp(['Calculating PLV for ' mat2str(sum(selectArr, 1)) ' trials...']);
for ichA = 1:nCh
    eegPhase(ichA, :, :) = angle(hilbert(squeeze(eegfilt(ichA, :, :))));
end

plv = zeros(size(eegPhase, 2), nCh, nCh, nCond);
for ichA = 1:nCh-1
    phaseA = squeeze(eegPhase(ichA, :, :));
    for ichB = ichA+1:nCh
        phaseB = squeeze(eegPhase(ichB, :, :));
        for iCond = 1:nCond
            plv(:, ichA, ichB, iCond) = abs(sum(exp(1i*(phaseA(:, selectArr(:, iCond)) ...
                                      - phaseB(:, selectArr(:, iCond)))), 2))/sum(selectArr(:, iCond));
        end
    end
end
plv = squeeze(plv);
return;


















