% Provide Eyelink with details about the graphics environment
% and perform some initializations. The information is returned
% in a structure that also contains useful defaults
% and control codes (e.g. tracker state bit and Eyelink key values).

    el = EyelinkInitDefaults(window);
    
    el.backgroundcolour        = BlackIndex(log.window);
    el.msgfontcolour           = WhiteIndex(log.window);
    el.imgtitlecolour          = WhiteIndex(log.window);
    el.calibrationtargetcolour = WhiteIndex(log.window);
    
    EyelinkUpdateDefaults(el);
    
    
% Before the first calibration


% Initialization of the connection with the Eyelink Gazetracker.
% exit program if this fails.
    if ~EyelinkInit(dummymode, 1)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end

% open file to record data to
%    edfFile=log.sonaNo;
    edfFile = strcat(log.sonaNo,'.edf');
    res = Eyelink('Openfile', edfFile);
    if res~=0
        fprintf('Cannot create EDF file ''%s'' ', edffilename);
        cleanup;
        return;
    end
    
% make sure we're still connected.
    if Eyelink('IsConnected')~=1 && ~dummymode
        cleanup;
        return;
    end
    
% This command is crucial to map the gaze positions from the tracker to
% screen pixel positions to determine fixation
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, winWidth-1, winHeight-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, winWidth-1, winHeight-1);
    
    
    [v, vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    vsn = regexp(vs,'\d','match');
    
    if v ==3 && str2double(vsn{1}) == 4 % if EL 1000 and tracker version 4.xx
% remote mode possible add HTARGET ( head target)
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT,HTARGET');
% set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT,HTARGET');
    else
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT');
% set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    end
    
% Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
% do a final check of calibration using driftcorrection
    success=EyelinkDoDriftCorrection(el);
    if success~=1
        cleanup;
        return;
    end
    
    
    
% The following should be pasted right at the start of a trial; perhaps
% at the top of the trial loop (if you use that :P)


% Sending a 'TRIALID' message to mark the start of a trial in Data 
% Viewer.  This is different than the start of recording message 
% START that is logged when the trial recording begins. The viewer
% will not parse any messages, events, or samples, that exist in 
% the data file prior to this message. 
    Eyelink('Message', 'TRIALID %d', t);
    Eyelink('Message', 'BLOCKID %d', block);

    Eyelink('command', 'record_status_message "TRIAL %d/%d, BLOCK %d/%d"', t, n, block, b);
            
% start recording eye position (preceded by a short pause so that 
% the tracker can finish the mode transition)
% The paramerters for the 'StartRecording' call controls the
% file_samples, file_events, link_samples, link_events availability
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.05);
% Eyelink('StartRecording', 1, 1, 1, 1);    
    Eyelink('StartRecording');
            
% record a few samples before we actually start displaying
% otherwise you may lose a few msec of data 
    WaitSecs(0.1);
            
            
% When presenting the different objects during a trial, you want to make
% sure you send a message to the eyetracker (it will put it in your output
% file), so that you know when the different trial components were
% presented when looking at your eyetracking data. Below is an example of
% such a message:

    Eyelink('Message', 'PROBE');



% At the very end of a trial (so right before the 'end' of your trial loop... again if you use that :P) you want to
% out the following lines of code:

% adds 100 msec of data to catch final events
    WaitSecs(0.1);
% stop the recording of eye-movements for the current trial
    Eyelink('StopRecording');
        
% IMPORTANT! Don't send too many messages in a very short period of
% time or the EyeLink tracker may not be able to write them all 
% to the EDF file.
% Consider adding a short delay every few messages.

    for m=1:set(setInd(trialInd(t)))
        Eyelink('Message', '!V IAREA ELLIPSE %d %d %d %d %d %s', m, xBCen(posInd(m))-lRad, yBCen(posInd(m))-lRad, xBCen(posInd(m))+lRad, yBCen(posInd(m))+lRad,'center');
        Eyelink('Message', '!V IAREA ELLIPSE %d %d %d %d %d %s', m, xBCen(posInd(m))-lRad*cos(randTheta(m)+pi/2)-pRad,yBCen(posInd(m))-lRad*sin(randTheta(m)+pi/2)-pRad,xBCen(posInd(m))-lRad*cos(randTheta(m)+pi/2)+pRad,yBCen(posInd(m))-lRad*sin(randTheta(m)+pi/2)+pRad,'center');
        WaitSecs(0.001);
    end
        
    Eyelink('Message', 'CONDITION %s', stimCond);
    Eyelink('Message', 'SETSIZE %d', set(setInd(trialInd(t))));
      
    WaitSecs(0.001);
% Sending a 'TRIAL_RESULT' message to mark the end of a trial in 
% Data Viewer. This is different than the end of recording message 
% END that is logged when the trial recording ends. The viewer will
% not parse any messages, events, or samples that exist in the data 
% file after this message.
    Eyelink('Message', 'TRIAL_RESULT 0');
    
%%% END TRIAL


% After the end of a block (so after the 'end' that shuts down the block loop, yet again: if you use that :P) you want
% to paste the following:

%%% END BLOCK

% End of Experiment; close the file first   
% close data file and shut down tracker
        
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
 

% Next you want to copy the EDF file from the Eyelink computer to the
% computer on which you ran the task with the following lines:

    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        status=Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
    catch %#ok<*CTCH>
        fprintf('Problem receiving data file ''%s''\n', edfFile );
    end

 
    
% At the very end of your file, paste the following line which shuts down
% the eyelink connection:

    Eyelink('Shutdown');
    

% That should be all :-)!