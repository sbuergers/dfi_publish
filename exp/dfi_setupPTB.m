function setup = dfi_setupPTB(setup, debug)
% Nested function to set up screen in dfi_present_stimuli
% only added as real function to test functionality.
%
% ----------------------------------------
% adapted from Agoston Mihalic, 
% last updated, April 2015
%
try
    % Is PTB using OpenGL?
    AssertOpenGL;
    
    stim = setup.stim;

    % --- Sound ---
    stim.aud.lat.lev  = 1; % really need low latency mode
    stim.aud.lat.bias = 0; %0.01106; % adjust for fixed hardware bias
    % (not used at the moment and needs to be calibrated with e.g. micro & photo diode)
    InitializePsychSound(stim.aud.lat.lev); % PsychPortAudio sound driver

    % Open default sound device
    % Office: Audio subsystem is Windows DirectSound, Audio device is Primary Sound Driver, device id is 8
    setup.haudio = PsychPortAudio('Open', [], [], [], stim.aud.freq, 2);

    % Control over volume and latency bias
    PsychPortAudio('Volume', setup.haudio, stim.aud.vol); 
    PsychPortAudio('LatencyBias', setup.haudio, stim.aud.lat.bias);
    
    setup.stim = stim;

    % --- Screen ---
    disp = setup.disp;
    if debug
        PsychDebugWindowConfiguration % use transparent screen
        Screen('Preference', 'SkipSyncTests', 1); % skip screen synchronization, otherwise error due to PTB debug mode
        disp.old.verb = Screen('Preference', 'Verbosity', 1); % print only errors
    else
        Screen('Preference', 'SkipSyncTests', 0); % synchronize screen
        disp.old.verb   = Screen('Preference', 'Verbosity'); % print error+warning+status messages
        disp.old.visdeb = Screen('Preference', 'VisualDebugLevel', 3); % turn off PsychToolbox Welcome Sign
        % hide mouse cursor for experiment
        if strcmp(setup.name, 'PiaCRT') || strcmp(setup.name, 'SteveCRT') && ~strcmp(setup.subj.id, '999')
            HideCursor
        end 
    end
    disp.h = max(Screen('Screens')); % Get screen id

    % Get second screen id and give PTB control over it
    % to avoid conflict with other software using the GPU
    if numel(Screen('Screens')) > 1
        disp.hmacbook= min(Screen('Screens')); 
        if ~debug, Screen('OpenWindow', disp.hmacbook); end
    end
    
    % colors
    [disp.white, disp.black] = deal(WhiteIndex(disp.h), BlackIndex(disp.h));
    disp.red   = [255 0 0];
    disp.green = [0 255 0];
    disp.grey  = disp.white * stim.vis.bgcontr;
    if disp.incolor 
        [disp.win,  disp.rect] = Screen('OpenWindow', disp.h, disp.green);
    else
        [disp.win,  disp.rect] = Screen('OpenWindow', disp.h, disp.grey);
    end
    [disp.xcen, disp.ycen] = RectCenter(disp.rect);
    disp.oldPriority = Priority(MaxPriority(disp.win));     % High priority
    disp.ifi         = Screen('GetFlipInterval', disp.win); % 1/refresh rate
    fprintf(['\nPTB estimates a refresh rate of:', num2str(disp.ifi), '\n\n'])
    % Set text font and size
    Screen(disp.win, 'TextFont', setup.text.font);
    Screen(disp.win, 'TextSize', setup.text.size);
    
    setup.disp = disp;
    
catch ME % This section is executed only in case an error happens
    sca
    ListenChar(1)
    warning('Could not set up Psychtoolbox screen!');
    rethrow(ME)
end 

return


