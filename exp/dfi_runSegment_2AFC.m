function [data, abort] = dfi_runSegment_2AFC(s, data, b)
%[data] = dfi_runSegment_2AFC(s, data, b)
% starts the sequence of stimuli of one block of the sifi paradigm (see
% dfi_run_experiment for more detailed information on the paradigm). This
% function contains the trial-loop.
%
% INPUT:
%        s      =   experimental structure (including all specs)
%        data   =   trial dataset (so far without response data)
%        b      =   current block number
%
% OUTPUT:
%        data   =   response data and trial descriptives
%
% -----
% adapted from Agoston Mihalik, 
% last updated, April 2015
%

    %
    %       Stimulus overview
    %
    %     STIM1   Stim2      ID
    %       V1    V1V2   =   2            V 
    %      V1V2    V1    =   3            V 
    %       A1    A1A2   =   4            A only
    %      V1A1   V2A1   =   5            Fus
    %      V2A1   V1A1   =   6            Fus 
    %      A1A2    A1    =   7            A 
    %      V1A2   V2A2   =   8            Fis
    %      V2A2   V1A2   =   9            Fis 
    %   
    
    
    
%++++++++++++++++++++
%% 1.)   SETUP
%++++++++++++++++++++

% return true if abort key is pressed
abort = 0;

% initialize variables

fix             = s.stim.fix;
sAudio          = zeros(s.ntrls,4);
trlid           = data.trlid;
trgid           = trlid;
disp            = s.disp;
w               = disp.win;

% generate ITIs (mean +- std)
plusminus = Shuffle([ones(s.ntrls,1); -ones(s.ntrls,1)]);
plusminus = plusminus(1:numel(plusminus)/2);
ITI       = s.ITI(1) + plusminus .* rand(s.ntrls, 1) * s.ITI(2);

% segment soas
soa     = data.soa(1:s.ntrls);  
soa2AFC = s.stim.SOA2AFC;

% for all intervals get exact multiples of ifi
ifi     = disp.ifi;
slack   = round(s.slack/ifi)     * ifi;
ITI     = round(ITI/ifi)         * ifi;
resw    = round(s.reswin/ifi)    * ifi;
soa     = round(soa/ifi)         * ifi;
soa2AFC = round(soa2AFC/ifi)     * ifi;
vdur    = round(s.stim.vdur/ifi) * ifi;

% buffer sound
PsychPortAudio('FillBuffer', s.haudio, s.stim.aud.loc);

% play a silent beep to avoid timing issues with the first sound stimulus
oldVolume = PsychPortAudio('Volume', s.haudio, 0);
PsychPortAudio('Start', s.haudio, 1);  
PsychPortAudio('Stop', s.haudio, 1);
PsychPortAudio('Volume', s.haudio, oldVolume);

% start recording key presses
KbQueueStart(s.kb.subj);

% keep track of saccades (no central fixation)
saccCount = 0;

% send trigger for start of block
dfi_sendTrigger(s, s.trig.sBlock)



%++++++++++++++++++++++++
%% 2.)   TRIAL LOOP
%++++++++++++++++++++++++

% alternate between straight and rotated fixcross from trial to trial
rotdeg = repmat([0 45], [1 s.ntrls/2])';
for i=1:s.ntrls
    
    %.....................
    %% 2a.) Eyetracker
    %.....................
    if s.el.online
        Eyelink('Message', 'TRIALID %d', i);     % send 'TRIALID'
        % This supplies the title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message "TRIAL %d/%d, BLOCK %d/%d"', i, s.ntrls, b, s.blocks.N);
    end
    
    %..................................
    %% 2b.) Stimulus presentation
    %..................................
    % Fixation cross (idle)
    Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));    
    if i==1
        vbl = Screen(w, 'Flip');
    else
        vbl = Screen(w, 'Flip', stimStopTime+resw-ifi/2);
    end

    % Fixation cross (ready) - measure alpha here 
    Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
    vbl = Screen(w, 'Flip', vbl+slack-ifi/2);
    dfi_sendTrigger(s, s.trig.siti + trgid(i)) 
    Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));   

    % Trial types    
    % A
    if ismember(trlid(i), [4,7])                                    % S1A       
        trltypeA{1} = 'V0A1';                    
    elseif ismember(trlid(i), [2,3])             
        trltypeA{1} = 'V1A0';                    
    elseif ismember(trlid(i), [5,6,8,9])        
        trltypeA{1} = 'V1A1';    
    else
        trltypeA{1} = 'Nil';     
    end                                        

    if ismember(trlid(i), [7,8])                                    % S2A            
        trltypeA{2} = 'A2';                      
    elseif ismember(trlid(i), [3,6])
        trltypeA{2} = 'V2';
    elseif ismember(trlid(i), 9)
        trltypeA{2} = 'V2A2';
    else
        trltypeA{2} = 'Nil';
    end
    % B                                       
    if ismember(trlid(i), [4,7])                                    % S1B             
        trltypeB{1} = 'V0A1';                    
    elseif ismember(trlid(i), [2,3])             
        trltypeB{1} = 'V1A0';                    
    elseif ismember(trlid(i), [5,6,8,9])        
        trltypeB{1} = 'V1A1'; 
    else
        trltypeB{1} = 'Nil';
    end                                        

    if ismember(trlid(i), [4,9])                                    % S2B            
        trltypeB{2} = 'A2';                      
    elseif ismember(trlid(i), [2,5])
        trltypeB{2} = 'V2';
    elseif ismember(trlid(i), 8)
        trltypeB{2} = 'V2A2';
    else
        trltypeB{2} = 'Nil';
    end


    % ---------------------------------------------------------------------
    %%                             STIMULI A
    % ---------------------------------------------------------------------
    switch trltypeA{1} 
        case 'V1A1'
            PsychPortAudio('Start', s.haudio, 1, vbl+ITI(i)+ifi/2+s.stim.vposdelay);                      % A1
            if ~strcmp(s.stim.vis.type, 'circle')
                Screen(w, 'DrawTexture', s.stim.vis.texture, [], s.stim.vis.rect);                          
            else
                Screen(w, 'FillOval', disp.white, s.stim.vis.rect);                                     
            end
            Screen(w, 'DrawingFinished');
            vbl = Screen(w, 'Flip', vbl+ITI(i)-ifi/2);                                                    % V1
            dfi_sendTrigger(s, s.trig.a.sAV + trgid(i))                                                   % Tr.AV
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
            Screen(w, 'DrawingFinished');
            Screen(w, 'Flip', vbl+vdur-ifi/2);                                                            % V1off
            sAudio(i,1) = PsychPortAudio('Stop', s.haudio, 1);
        case 'V1A0' 
            if ~strcmp(s.stim.vis.type, 'circle')
                Screen(w, 'DrawTexture', s.stim.vis.texture, [], s.stim.vis.rect);                        % V1buff
            else
                Screen(w, 'FillOval', disp.white, s.stim.vis.rect);                                     
            end                               
            Screen(w, 'DrawingFinished');
            vbl = Screen(w, 'Flip', vbl+ITI(i)-ifi/2);                                                    % V1
            dfi_sendTrigger(s, s.trig.a.sV + trgid(i))                                                    % Tr.V
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
            Screen(w, 'DrawingFinished');
            Screen(w, 'Flip', vbl+vdur-ifi/2);                                                            % V1off
        case 'V0A1'
            PsychPortAudio('Start', s.haudio, 1, vbl+ITI(i)+ifi/2+s.stim.vposdelay);                      % A1
            dfi_sendTrigger(s, s.trig.a.sA + trgid(i))                                                    % Tr.A
            sAudio(i,1) = PsychPortAudio('Stop', s.haudio, 1);
        case 'Nil'
            disp('Stimulus:: None')
            dfi_sendTrigger(s, s.trig.a.nil + trgid(i))                                                   % Tr.noStim
            WaitSecs(vdur); 
    end

    % Display V2/A2
    switch trltypeA{2} 
        case 'V2A2' 
            PsychPortAudio('Start', s.haudio,1,sAudio(i,1)+soa(i));                                       % A2
            if ~strcmp(s.stim.vis.type,  'circle')
                Screen(w, 'DrawTextures', [s.stim.vis.texture, fix.texture], [], ...                      % V2buff
                                          [s.stim.vis.rect;    fix.rect]', [0, rotdeg(i)]);
            else
                Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i)); 
                Screen(w, 'FillOval', disp.white, s.stim.vis.rect);                                     
            end
            Screen(w, 'DrawingFinished');
            vbl = Screen(w, 'Flip', vbl+soa(i)-ifi/2);                                                    % V2
            dfi_sendTrigger(s, s.trig.a.sAV2 + trgid(i))                                                  % Tr.AV2
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
            Screen(w, 'DrawingFinished');
            Screen(w, 'Flip', vbl+vdur-ifi/2);                                                            % V2off
            sAudio(i,2) = PsychPortAudio('Stop', s.haudio, 1);
        case 'V2'
            if ~strcmp(s.stim.vis.type,  'circle')
                Screen(w, 'DrawTextures', [s.stim.vis.texture, fix.texture], [], ...                      % V2buff
                                          [s.stim.vis.rect;    fix.rect]', [0,rotdeg(i)]);    
            else
                Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i)); 
                Screen(w, 'FillOval', disp.white, s.stim.vis.rect);                                     
            end                                        
            Screen(w, 'DrawingFinished');
            vbl = Screen(w, 'Flip', vbl+soa(i)-ifi/2);                                                    % V2
            dfi_sendTrigger(s, s.trig.a.sV2 + trgid(i))                                                   % Tr.V2
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
            Screen(w, 'DrawingFinished');
            Screen(w, 'Flip', vbl+vdur-ifi/2);                                                            % V2off
        case 'A2'
            PsychPortAudio('Start', s.haudio,1,sAudio(i,1)+soa(i));                                       % A2
            dfi_sendTrigger(s, s.trig.a.sA2 + trgid(i))                                                   % Tr.A2
            sAudio(i,2) = PsychPortAudio('Stop', s.haudio, 1);
        otherwise
            WaitSecs(soa(i)+vdur); 
    end
    

    % ---------------------------------------------------------------------
    %%                              STIMULI B
    % ---------------------------------------------------------------------  
    Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
    vbl = Screen(w, 'Flip');  
    Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
    switch trltypeB{1} 
        case 'V1A1'
            PsychPortAudio('Start', s.haudio, 1, vbl+soa2AFC+ifi/2+s.stim.vposdelay);                     % A1
            if ~strcmp(s.stim.vis.type, 'circle')
                Screen(w, 'DrawTexture', s.stim.vis.texture, [], s.stim.vis.rect);                          
            else
                Screen(w, 'FillOval', disp.white, s.stim.vis.rect);                                     
            end
            Screen(w, 'DrawingFinished');
            vbl = Screen(w, 'Flip', vbl+soa2AFC-ifi/2);                                                   % V1
            dfi_sendTrigger(s, s.trig.b.sAV + trgid(i))                                                   % Tr.AV
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
            Screen(w, 'DrawingFinished');
            Screen(w, 'Flip', vbl+vdur-ifi/2);                                                            % V1off
        case 'V1A0' 
            if ~strcmp(s.stim.vis.type, 'circle')
                Screen(w, 'DrawTexture', s.stim.vis.texture, [], s.stim.vis.rect);                        % V1buff
            else
                Screen(w, 'FillOval', disp.white, s.stim.vis.rect);                                     
            end                               
            Screen(w, 'DrawingFinished');
            vbl = Screen(w, 'Flip', vbl+soa2AFC-ifi/2);                                                    % V1
            dfi_sendTrigger(s, s.trig.b.sV + trgid(i))                                                    % Tr.V
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
            Screen(w, 'DrawingFinished');
            Screen(w, 'Flip', vbl+vdur-ifi/2);                                                            % V1off
        case 'V0A1'
            PsychPortAudio('Start', s.haudio, 1, vbl+soa2AFC+ifi/2+s.stim.vposdelay);                     % A1
            dfi_sendTrigger(s, s.trig.b.sA + trgid(i))                                                    % Tr.A
        case 'nil'
            disp('Stimulus:: None')
            dfi_sendTrigger(s, s.trig.b.nil + trgid(i))                                                   % Tr.noStim
            WaitSecs(vdur); 
    end
    sAudio(i,3) = PsychPortAudio('Stop', s.haudio, 1);                                                    % A1off
    % Display V2/A2
    switch trltypeB{2} 
        case 'V2A2' 
            PsychPortAudio('Start', s.haudio,1,sAudio(i,3)+soa(i));                                       % A2
            if ~strcmp(s.stim.vis.type,  'circle')
                Screen(w, 'DrawTextures', [s.stim.vis.texture, fix.texture], [], ...                      % V2buff
                                          [s.stim.vis.rect;    fix.rect]', [0, rotdeg(i)]);
            else
                Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i)); 
                Screen(w, 'FillOval', disp.white, s.stim.vis.rect);                                     
            end
            Screen(w, 'DrawingFinished');
            vbl = Screen(w, 'Flip', vbl+soa(i)-ifi/2);                                                    % V2
            dfi_sendTrigger(s, s.trig.b.sAV2 + trgid(i))                                                  % Tr.AV2
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
            Screen(w, 'DrawingFinished');
            Screen(w, 'Flip', vbl+vdur-ifi/2);                                                            % V2off
        case 'V2'
            if ~strcmp(s.stim.vis.type,  'circle')
                Screen(w, 'DrawTextures', [s.stim.vis.texture, fix.texture], [], ...                      % V2buff
                                          [s.stim.vis.rect;    fix.rect]', [0, rotdeg(i)]);                
            else
                Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i)); 
                Screen(w, 'FillOval', disp.white, s.stim.vis.rect);                                     
            end                                        
            Screen(w, 'DrawingFinished');
            vbl = Screen(w, 'Flip', vbl+soa(i)-ifi/2);                                                    % V2
            dfi_sendTrigger(s, s.trig.b.sV2 + trgid(i))                                                   % Tr.V2
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
            Screen(w, 'DrawingFinished');
            Screen(w, 'Flip', vbl+vdur-ifi/2);                                                            % V2off
        case 'A2'
            PsychPortAudio('Start', s.haudio,1,sAudio(i,3)+soa(i));                                       % A2
            dfi_sendTrigger(s, s.trig.b.sA2 + trgid(i))                                                   % Tr.A2
        otherwise
            WaitSecs(soa(i)); 
    end
    sAudio(i,4) = PsychPortAudio('Stop', s.haudio, 1);                                                    % A2off

    % ************************* END OF STIMULI ************************** %

    
    %...............................
    %% 2c.) Response collection
    %...............................

    % Have there been responses in the meantime?
    [event_mult, nremaining] = KbEventGet(s.kb.subj);
    if nremaining > 0 && i > 1
        if exist('stimStopTime', 'var')
            data.badRT(i-1)    = event_mult.Time - stimStopTime;
        else
            data.badRT(i-1)    = -99;
        end
        data.nbadresp(i-1) = nremaining/2;
        if event_mult.Keycode == KbName(s.key.quit)
            abort = 1;
            if s.el.online
                Eyelink('Message', 'RUN ABORTED BY USER');
            end
            return
        end
    end

    % resynchronize to vbl
    Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
    stimStopTime = Screen(w, 'Flip');

    % save response time also relative to current trial
    if nremaining > 0
        data.badRT(i)    = event_mult.Time - stimStopTime;
        data.nbadresp(i) = KbEventFlush(s.kb.subj)/2; 
    end

    switch strcmp(s.respdev, 'keyboard')
        
        % collect responses with the keyboard
        case 1 % keyboard
            while data.resp(i)==0 && GetSecs < stimStopTime + resw-ifi
                [event, nremaining] = KbEventGet(s.kb.subj);
                % there is a nice single response
                if isstruct(event) && event.Pressed == 1 && nremaining <= 1
                    if ismember(event.Keycode, s.key.code)
                        data.RT(i)   = event.Time - stimStopTime;
                        data.resp(i) = find(event.Keycode==s.key.code);
                        % send EEG & EL triggers
                        if data.resp(i) == 1
                            dfi_sendTrigger(s, s.trig.rOne + trgid(i))
                        else
                            dfi_sendTrigger(s, s.trig.rTwo + trgid(i))
                        end
                    end
                    KbEventFlush(s.kb.subj)
                end
            end % end response collection
            
        % collect responses with the mouse (includes confidence judgment)
        case 0 % mouse
            % show cursor in center
            SetMouse(disp.xcen, disp.ycen, w); 
            ShowCursor;
            
            % draw butterfly confidence judgment window
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
            ab = dfi_conf_butterfly_draw(s);
            Screen(w, 'Flip');
            
% stimStopTime = GetSecs % uncomment for debugging...
            while data.resp(i)==0 && GetSecs < stimStopTime + resw-ifi
                
                % get mouse coordinates and clicks
                button_status = 0; % cursor not in butterfly
                [x, y, buttons] = GetMouse(w);
                tclick = GetSecs;
                ab.x = x; 
                ab.y = y;
                
                % draw fixation cross and confidence butterfly
                Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
                dfi_conf_butterfly_draw(s);
                
                % highlight buttons when mouse is hovering over them, 
                % find angle relative to vertical axis
                opp = y - disp.ycen;
                hyp = sqrt(opp.^2 + (x - disp.xcen).^2);
                ang = -asin(opp/hyp)/pi*180; % negative for low screen
                % is the angle in the butterfly range (in degrees)?
                if abs(ang) < ab.b_angle*ab.b_num/2
                    % get hypothenuse in cm not pixels
                    ppcmw  = disp.res(1)/disp.width;
                    ppcmh  = disp.res(2)/disp.height;
                    hyp_cm = sqrt(  (abs(y-disp.ycen) / ppcmh).^2 + (abs(x-disp.xcen) / ppcmw).^2  );
                    % get hypothenuse in visual angle
                    % in this case the hypothenuse is actually the opposite...
                    hyp_ang = atan(hyp_cm / disp.dist) * 180/pi;
                    % is the response within the response butterfly?
                    if ab.inner_diam/2 < hyp_ang && hyp_ang < ab.outer_diam/2
                        % fill response button
                        ab.resp_angle = ang;
                        button_status = dfi_conf_butterfly_fill(s, ab);
                        Screen('DrawTexture', disp.win, fix.texture, [], fix.rect, rotdeg(i));
                    end
                end
                % Flip the Screen
                Screen('DrawingFinished', w);
                Screen(w, 'Flip');
                
                % save data if a button was clicked
                if any(buttons) && button_status ~= 0
                    data.RT(i)   = tclick - vbl+ifi/2;
                    if button_status > 0
                        data.resp(i) = 2;
                    elseif button_status < 0
                        data.resp(i) = 1;
                    end
                    data.conf(i) = abs(button_status);
%                     Screen(w, 'DrawTexture', fix.texture, [], fix.rect, rotdeg(i));
%                     Screen(w, 'Flip');
                end
            end
            HideCursor;
    end % switch respdev

    % mark end of trial in eyelink log
    if s.el.online
        Eyelink('Message', 'END TRIAL');
        [saccadeOccured] = dfi_checkForSaccades(s);
        if saccadeOccured
            saccCount = saccCount + 1;
            % Feedback about central fixation after a trial
            if strcmp(s.eyelink.fb, 'trial')
                Screen(disp.win, 'Flip'); 
                DrawFormattedText(disp.win, 'Fixate on the center', 'center', 'center', disp.white);
                Screen(disp.win, 'Flip'); 
                WaitSecs(1.5); 
                Screen(disp.win, 'DrawTexture', s.stim.fix.texture, [], s.stim.fix.rect, rotdeg(i));
                Screen(disp.win, 'Flip'); 
            end
        end
    end
    
    % pre-practice feedback
    if s.practicefeedback 
        if strcmp(data.rkeycfg(1), 'AB')
            if data.resp(i)==2 && ismember(data.trlid(i), [2,4,5,8]) || ...
               data.resp(i)==1 && ismember(data.trlid(i), [3,6,7,9])
                dfi_instructions('exampleFeedbackCorr', s)
            else
                dfi_instructions('exampleFeedbackIncorr', s)
            end
        elseif strcmp(data.rkeycfg(1), 'BA')
            if data.resp(i)==1 && ismember(data.trlid(i), [2,4,5,8]) || ...
               data.resp(i)==2 && ismember(data.trlid(i), [3,6,7,9])
                dfi_instructions('exampleFeedbackCorr', s)
            else
                dfi_instructions('exampleFeedbackIncorr', s)
            end
        end
    end
    
end % end trial loop 

% Feedback about central fixation after a block
if saccCount > 1 && strcmp(s.eyelink.fb, 'block')
    Screen(disp.win, 'Flip'); 
    DrawFormattedText(disp.win, 'Please try to fixate on the cross in the middle of the screen.', 'center', 'center', disp.white);
    Screen(disp.win, 'Flip'); 
    WaitSecs(3);  
end


%+++++++++++++++++++++
%% 3.)   Finish
%+++++++++++++++++++++

% send trigger for end of block
dfi_sendTrigger(s, s.trig.eBlock)

% Empty keyboard queue
KbQueueFlush(s.kb.subj);

% no vbl times here retained. For this check dfi_runSegment, which is the
% yesno version of this task. The code is basically the same and there is
% no reason to assume that the flip diagnostics should be different between
% these two versions except for the interval between the two stimulus sets,
% and there it doesn't matter very much if a flip is missed or not.

% sAudio - sAudio(1,1)
% a = ...
%         [0    0.0267    0.2405    0.2672; ...
%     1.0945    1.1211    1.3617    1.3884; ...
%     2.1890    2.2690    2.4829    2.5629; ...
%     3.3903    3.6438    3.8577    3.8577; ...
%     4.9386    4.9386    5.2326    5.2326; ...
%     6.1134    6.1667    6.3806    6.4339; ...
%     7.2612    7.3145    7.5285    7.5819; ...
%     8.4092    8.4626    8.6764    8.6764; ...
%     9.5570    9.6638    9.9043   10.0110; ...
%    10.8117   10.9184   11.1591   11.1591; ...
%    12.0931   12.1731   12.4137   12.4137; ...
%    13.3211   13.3745   13.6150   13.6684; ...
%    14.4690   14.4690   14.7362   14.7629; ...
%    15.5636   15.5636   15.8575   15.8575; ...
%    16.7382   16.7382   17.0055   17.0055; ...
%    17.8594   17.8861   18.1266   18.1533; ...
%    18.9540   19.0607   19.2745   19.3812; ...
%    20.2086   20.2086   20.5292   20.6093; ...
%    21.4366   21.4366   21.7305   21.7839; ...
%    22.6112   22.8647   23.0787   23.3321; ...
%    24.1596   24.1863   24.4002   24.4002; ...
%    25.2541   25.2541   25.7482   25.7482; ...
%    26.8292   27.0827   27.3234   27.5769; ...
%    28.3776   28.3776   28.6448   28.6715; ...
%    29.4721   29.4721   29.8194   29.9262; ...
%    30.7268   30.8068   31.0474   31.1274; ...
%    31.9281   31.9548   32.1687   32.1953; ...
%    33.0227   33.0760   33.3165   33.3699];
%++++++++++++++++++++++++++++
% 4.)   Nested functions
%++++++++++++++++++++++++++++

function dfi_sendTrigger( s, trig )
%sendTrigger 
% uses a LabJack device to send triggers to an external hardware, like an
% EEG setup. Needs labjack object (s.lj) and the trigger value to be send (trig). 
    % Send trigger now. 
    if s.ljPresent
        prepareStrobe(s.lj, trig);
        strobeWord(s.lj);
    end
    % Send eyelink trigger
    if s.el.online
        Eyelink('Message', ['TRIGGER ', num2str(trig)]);
    end
end


% unfinished nested functions...
% ...


% function dfi_checkGazePosition( s )
% % did the gaze position change? If so tell me how!
%     if Eyelink( 'NewFloatSampleAvailable') 
%         % get the sample in the form of an event structure
%         evt = Eyelink( 'NewestFloatSample');
%         if s.el.eye_used ~= -1 % do we know which eye to use yet?
%             % if we do, get current gaze position from sample
%             x = evt.gx(s.el.eye_used+1); % +1 as we're accessing MATLAB array
%             y = evt.gy(s.el.eye_used+1);
%             % do we have valid data and is the pupil visible?
%             if x~=s.el.MISSING_DATA && y~=s.el.MISSING_DATA && evt.pa(s.el.eye_used+1)>0
%                 s.el.mx=x;
%                 s.el.my=y;
%             end
%         end
%     end
% end

% 
% % check for endsaccade events
% if Eyelink('isconnected') == el.dummyconnected % in dummy mode use mousecoordinates
%     [x,y,button] = GetMouse(window);
%     evt.type=el.ENDSACC;
%     evt.genx=x;
%     evt.geny=y;
%     evtype=el.ENDSACC;
% else % check for events
%     evtype=Eyelink('getnextdatatype');
% end
% if evtype==el.ENDSACC		% if the subject finished a saccade check if it fell on an object
%     if Eyelink('isconnected') == el.connected % if we're really measuring eye-movements
%         evt = Eyelink('getfloatdata', evtype); % get data
%     end
%     % check if saccade landed on an object
%     choice=-1;
%     noobject=0;
%     i=1;
%     while 1
%         if 1==IsInRect(evt.genx,evt.geny, object(i).rect )
%             choice=i;
%             break;
%         end
%         i=i+1;
%         if i>length(object)
%             noobject=1;
%             break;
%         end
%     end
%     if lastchoice>0 && (choice~=lastchoice || noobject==1) % toggle object color
%         if object(lastchoice).on==1 % restore screen
%             Screen('CopyWindow', buffer, window, object(lastchoice).rect, object(lastchoice).rect);
%             object(lastchoice).on=0;
%             lastchoice=-1;
%             doflip=1;
%         end
%     end
%     if choice>0 && choice~=lastchoice % toggle object color
%         if object(choice).on==0 % toggle object on screen
%             Screen('CopyWindow', altbuffer, window, object(choice).rect, object(choice).rect);
%             object(choice).on=1;
%             doflip=1;
%         end
%         lastchoice=choice;
%     end
%     if doflip==1
%         Screen('Flip',  window, [], 1);
%         doflip=0;
%     end
% end % saccade?


end % eof
    
    
    
    
    
    
    
    
    
    
    