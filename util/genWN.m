function wn = genWN(Fs, stim)
    
    % gaussian white noise in audiosamples of stimulus duration
    wncen = randn(round(Fs * stim.adur), 1) / 3; 
    % copy for both ears
    wn = zeros(Fs*stim.adur, 2, 1);
    wn(:,1,:) = wncen;
    wn(:,2,:) = wncen;
    wn(wn> 1) =  1; % cutoff above +3 std
    wn(wn<-1) = -1; % cutoff below -3 std 

end % end wn