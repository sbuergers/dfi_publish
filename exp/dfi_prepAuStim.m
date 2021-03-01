function stim = dfi_prepAuStim(stim)
% auditory pure tone stimulus with ramp on-off for fast, 
% accurate presentation.
%
% ----------------------------------------
% adapted from
% http://www.h6.dion.ne.jp/~fff/old/technique/auditory/matlab.html
%    

    % prepare tone
    cf = stim.afreq;             % carrier frequency (Hz)
    sf = stim.aud.freq;          % sample frequency (Hz)
    Fs = sf;
    d  = stim.adur;              % duration (s)
    tlag = stim.ITD.timelag; % ITD in seconds
    slag = round(sf * tlag); % ITD in samples by PortAudio
     
    % use a white noise stimulus to apply HRTF to
    aloc = genWN(Fs, stim);
        
    
    % adjust volume
    aloc = aloc * stim.aud.vol;
    
    if stim.plot
        subplot(2,1,1)
        sound(aloc, sf);            % sound presentation
        pause(d + 0.5);             % waiting for sound end
        plot(aloc);
        xlabel('sound samples')     % x-axis label
        ylabel('sine values')       % y-axis label
    end
        
    % note that I cannot use the HRTF because the convolution 
    % requires a window of a certain size, which I do not have.
    % So use ITD instead
    if stim.vloc(1) < 0
        aloc(slag+1:end,1) = aloc(1:end-slag,1); % left
        aloc(1:slag,1)     = 0;
    elseif stim.vloc(1) > 0
        aloc(slag+1:end,2) = aloc(1:end-slag,2); % right
        aloc(1:slag,2)     = 0;
    end
    
    
    % ramp?
    if stim.aud.ramp
        r = round(Fs * stim.aud.rdur); % calculate r
        aloc(1:r,:) = aloc(1:r,:) .* repmat((1:r)',[1,2]) / r; % add r on
        aloc(end-r+1:end,:) = aloc(end-r+1:end,:) .* repmat((r:-1:1)', [1,2]) / r; % add r off
    end
    
    if stim.plot
        subplot(2,1,2)
        sound(aloc, sf);            % sound presentation
        pause(d + 0.5);             % waiting for sound end
        plot(aloc(:,1), 'g'); hold on
        plot(aloc(:,2), 'b');
        legend('left', 'right')
        xlabel('sound samples')     % x-axis label
        ylabel('sine values')       % y-axis label
    end
    
    stim.aud.loc = aloc'; % output
    
    
    %%  *** Nested functions ***
    
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


end % end dfi_prepAuStim


% eof








