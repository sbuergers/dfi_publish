function [time, data, abort, s] = dfi_runSegment(s, data, b)
%[time, data, abort, s] = dfi_runSegment(s, data, b)
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
%        time   =   matrix with flip timestamps
%       abort   =   boolean (abort key was pressed)
%           s   =   updated experiment structure 's' including staircase
%                   information
%
% -----
% adapted from Agoston Mihalik, 
% last updated, Aug. 2016 (staircase)
%


%+++++++++++++++++++
%% 1.)   SETUP
%+++++++++++++++++++

% return true if abort key is pressed
time  = []; % default value in case q is pressed
abort = 0;
ntrls = s.ntrls;

% initialize variables
w               = s.disp.win;
fix             = s.stim.fix;
sAudio          = zeros(ntrls,2);
trlid           = data.trlid;
trgid           = trlid;

% generate ITIs (mean +- std)
plusminus = Shuffle([ones(ntrls,1); -ones(ntrls,1)]);
plusminus = plusminus(1:numel(plusminus)/2);
ITI       = s.ITI(1) + plusminus .* rand(ntrls, 1) * s.ITI(2);

% segment soas
soa = data.soa(1:ntrls);  

% setup flip-times matrix (time x nflips)
% VBLTimestamp StimulusOnsetTime FlipTimestamp Missed (Beampos)
vbltm = zeros(ntrls, 6); 
sotm  = zeros(ntrls, 6);
ftsm  = zeros(ntrls, 6);
missm = zeros(ntrls, 6);

% for all intervals get exact multiples of ifi
ifi   = s.disp.ifi;
slack = round(s.slack/ifi)     * ifi;
ITI   = round(ITI/ifi)         * ifi;
resw  = round(s.reswin/ifi)    * ifi;
soa   = round(soa/ifi)         * ifi;
vdur  = round(s.stim.vdur/ifi) * ifi;

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


%+++++++++++++++++++++++
%% 2.)   TRIAL LOOP
%+++++++++++++++++++++++

for i=1:ntrls
    
    %.......................
    %% Staircase selection
    %.......................
    if strcmp(s.stair.incl, 'yes')
        % randomly sample the chainID for every trial (with replacement!)
        % but only when the chain hasn't converged yet
        converged_already = 1;
        cnt = 0;
        while converged_already && cnt <= s.stair.nChains
            converged_already = 0;
            cnt = cnt + 1;
            chainID = randsample(s.stair.nChains ,1);
            if ismember(trlid(i), [2,3])
                soa(i) = s.ud{chainID}.FF.xCurrent;  
                if s.ud{chainID}.FF.stop == 1,  converged_already = 1; end;
                sprintf('\n\n\nChain %i, condition FF, SOA = %f2\n', chainID, soa(i))
            elseif ismember(trlid(i), [5,6])
                soa(i) = s.ud{chainID}.Fus.xCurrent; 
                if s.ud{chainID}.Fus.stop == 1, converged_already = 1; end;
                sprintf('\n\n\nChain %i, condition Fusion, SOA = %f2\n', chainID, soa(i))
            elseif ismember(trlid(i), [8,9])
                soa(i) = s.ud{chainID}.Fis.xCurrent; 
                if s.ud{chainID}.Fis.stop == 1, converged_already = 1; end;
                sprintf('\n\n\nChain %i, condition Fission, SOA = %f2\n', chainID, soa(i))
            end;
        end
        data.soa(i) = soa(i);
        data.chainid(i) = chainID;
        % go to next trial if chains of this condition have already converged.
        if converged_already
            % resynchronize to vbl
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect);
            continue
        end
    end
    
    
    %....................
    %% 2a.) Eyetracker
    %....................
    if s.el.online
        Eyelink('Message', 'TRIALID %d', i);     % send 'TRIALID'
        % This supplies the title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message "TRIAL %d/%d, BLOCK %d/%d"', i, ntrls, b, s.blocks.N);
    end
    
    
    %..................................
    %% 2b.) Stimulus presentation
    %..................................
    % Fixation cross (idle)
    Screen(w, 'DrawTexture', fix.texture, [], fix.rect);
    if i==1 || ~exist('stimStopTime', 'var')
        [vbltm(i,1), sotm(i,1), ftsm(i,1)] = Screen(w, 'Flip');
    else
        [vbltm(i,1), sotm(i,1), ftsm(i,1)] = Screen(w, 'Flip', stimStopTime+resw-ifi/2);
    end

    % Fixation cross (ready) - measure alpha here 
    Screen(w, 'DrawTexture', fix.texture, [], fix.rect);
    [vbltm(i,2), sotm(i,2), ftsm(i,2)] = Screen(w, 'Flip', vbltm(i,1)+slack-ifi/2);
    dfi_sendTrigger(s, s.trig.siti + trgid(i)) 
    Screen(w, 'DrawTexture', fix.texture, [], fix.rect);

    % trial type?
    % S1                                        % Stimulus overview
    if ismember(trlid(i), [4,7])                %     V A
        trltype{1} = 'V0A1';                    % 1   0 0
    elseif ismember(trlid(i), [2,3])            % 2   1 0 
        trltype{1} = 'V1A0';                    % 3   2 0
    elseif ismember(trlid(i), [5,6,8,9])        % 4   0 1
        trltype{1} = 'V1A1';                    % 5   1 1
    else                                        % 6   2 1
        trltype{1} = 'Nil';                     % 7   0 2
    end                                         % 8   1 2
    % S2                                        % 9   2 2
    if ismember(trlid(i), [7,8])                
        trltype{2} = 'A2';                      
    elseif ismember(trlid(i), [3,6])
        trltype{2} = 'V2';
    elseif ismember(trlid(i), 9)
        trltype{2} = 'V2A2';
    else
        trltype{2} = 'Nil';
    end
    

    % Display V1/A1
    % Good onset asynchrony for V and A? Simultaneous? (yes, Cecere)
    switch trltype{1} 
        case 'V1A1'
            PsychPortAudio('Start', s.haudio, 1, vbltm(i,2)+ITI(i)+ifi/2+s.stim.vposdelay);               % A1
            if ~strcmp(s.stim.vis.type, 'circle')
                Screen(w, 'DrawTexture', s.stim.vis.texture, [], s.stim.vis.rect);                          
            else
                Screen(w, 'FillOval', s.disp.white, s.stim.vis.rect);                                     
            end
            Screen(w, 'DrawingFinished');
            [vbltm(i,3), sotm(i,3), ftsm(i,3)] = Screen(w, 'Flip', vbltm(i,2)+ITI(i)-ifi/2);              % V1
            dfi_sendTrigger(s, s.trig.sAV + trgid(i))                                                     % Tr.AV
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect);
            Screen(w, 'DrawingFinished');
            [vbltm(i,4), sotm(i,4), ftsm(i,4)] = Screen(w, 'Flip', vbltm(i,3)+vdur-ifi/2);                % V1off
        case 'V1A0' 
            if ~strcmp(s.stim.vis.type, 'circle')
                Screen(w, 'DrawTexture', s.stim.vis.texture, [], s.stim.vis.rect);                        % V1buff
            else
                Screen(w, 'FillOval', s.disp.white, s.stim.vis.rect);                                     
            end                               
            Screen(w, 'DrawingFinished');
            [vbltm(i,3), sotm(i,3), ftsm(i,3)] = Screen(w, 'Flip', vbltm(i,2)+ITI(i)-ifi/2);              % V1
            dfi_sendTrigger(s, s.trig.sV + trgid(i))                                                      % Tr.V
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect);
            Screen(w, 'DrawingFinished');
            [vbltm(i,4), sotm(i,4), ftsm(i,4)] = Screen(w, 'Flip', vbltm(i,3)+vdur-ifi/2);                % V1off
        case 'V0A1'
            PsychPortAudio('Start', s.haudio, 1, vbltm(i,2)+ITI(i)+ifi/2+s.stim.vposdelay);               % A1
            dfi_sendTrigger(s, s.trig.sA + trgid(i))                                                      % Tr.A
        case 'Nil'
            disp('Stimulus:: None')
            dfi_sendTrigger(s, s.trig.nil + trgid(i))                                                     % Tr.noStim
            WaitSecs(vdur); 
            
    end
    sAudio(i,1) = PsychPortAudio('Stop', s.haudio, 1);                                                    % A1off
    % Display V2/A2
    switch trltype{2} 
        case 'V2A2' 
            PsychPortAudio('Start', s.haudio,1,sAudio(i,1)+soa(i));                                       % A2
            if ~strcmp(s.stim.vis.type,  'circle')
                Screen(w, 'DrawTextures', [s.stim.vis.texture, fix.texture], [], ...                      % V2buff
                                          [s.stim.vis.rect;    fix.rect]');
            else
                Screen(w, 'DrawTexture', fix.texture, [], fix.rect); 
                Screen(w, 'FillOval', s.disp.white, s.stim.vis.rect);                                     
            end
            Screen(w, 'DrawingFinished');
            [vbltm(i,5), sotm(i,5), ftsm(i,5)] = Screen(w, 'Flip', vbltm(i,3)+soa(i)-ifi/2);              % V2
            dfi_sendTrigger(s, s.trig.sAV2 + trgid(i))                                                    % Tr.AV2
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect);
            Screen(w, 'DrawingFinished');
            [vbltm(i,6), sotm(i,6), ftsm(i,6)] = Screen(w, 'Flip', vbltm(i,5)+vdur-ifi/2);                % V2off
        case 'V2'
            if ~strcmp(s.stim.vis.type,  'circle')
                Screen(w, 'DrawTextures', [s.stim.vis.texture, fix.texture], [], ...                      % V2buff
                                          [s.stim.vis.rect;    fix.rect]');                
            else
                Screen(w, 'DrawTexture', fix.texture, [], fix.rect); 
                Screen(w, 'FillOval', s.disp.white, s.stim.vis.rect);                                     
            end                                        
            Screen(w, 'DrawingFinished');
            [vbltm(i,5), sotm(i,5), ftsm(i,5)] = Screen(w, 'Flip', vbltm(i,3)+soa(i)-ifi/2);              % V2
            dfi_sendTrigger(s, s.trig.sV2 + trgid(i))                                                     % Tr.V2
            Screen(w, 'DrawTexture', fix.texture, [], fix.rect);
            Screen(w, 'DrawingFinished');
            [vbltm(i,6), sotm(i,6), ftsm(i,6)] = Screen(w, 'Flip', vbltm(i,5)+vdur-ifi/2);                % V2off
        case 'A2'
            PsychPortAudio('Start', s.haudio,1,sAudio(i,1)+soa(i));                                       % A2
            dfi_sendTrigger(s, s.trig.sA2 + trgid(i))                                                     % Tr.A2
        case 'Nil'
            WaitSecs(vdur+soa(i));
        otherwise
            WaitSecs(ITI(i)+2*vdur); 
    end
    sAudio(i,2) = PsychPortAudio('Stop', s.haudio, 1);                                                    % A2off

    
    %...............................
    %% 2c.) Response collection
    %...............................
    % Have there been responses in the meantime?
    [event_mult, nremaining_mult] = KbEventGet(s.kb.subj);
    if nremaining_mult > 0 && i > 1
        if exist('stimStopTime', 'var')
            data.badRT(i-1)    = event_mult.Time - stimStopTime;
        else
            data.badRT(i-1)    = -99;
        end
        data.nbadresp(i-1) = nremaining_mult/2;
        if event_mult.Keycode == KbName(s.key.quit)
            abort = 1;
            if s.el.online
                Eyelink('Message', 'RUN ABORTED BY USER');
            end
            return
        end
    end
    
    % resynchronize to vbl
    Screen(w, 'DrawTexture', fix.texture, [], fix.rect);
    stimStopTime = Screen(w, 'Flip');

    % save response time also for current trial
    if nremaining_mult > 0
        data.badRT(i)    = event_mult.Time - stimStopTime;
        data.nbadresp(i) = KbEventFlush(s.kb.subj)/2; 
    end
    
    % draw butterfly confidence judgment window
    if s.practice && strcmp(s.paradigm, 'YN_threshold')
        Screen(w, 'DrawTexture', fix.texture, [], fix.rect);
        ab = dfi_conf_butterfly_draw(s);
        Screen(w, 'Flip');
    end
    
    while data.resp(i)==0 && GetSecs < stimStopTime + resw-ifi
        [event, nremaining] = KbEventGet(s.kb.subj);
        % there is a nice single response
        if isstruct(event) && event.Pressed == 1 && nremaining <= 1
            % abort key pressed?
            if event.Keycode == KbName(s.key.quit)
                abort = 1;
                if s.el.online
                    Eyelink('Message', 'RUN ABORTED BY USER');
                end
                return
            end
            if ismember(event.Keycode, s.key.code)
                data.RT(i)   = event.Time - stimStopTime;
                data.resp(i) = find(event.Keycode==s.key.code);
                % send EEG triggers
                if data.resp(i) == 1
                    dfi_sendTrigger(s, s.trig.rOne + trgid(i))
                else
                    dfi_sendTrigger(s, s.trig.rTwo + trgid(i))
                end
                % fill response button (only for the first practice block)
                if s.practice && strcmp(s.paradigm, 'YN_threshold') && b == 1
                    switch data.resp(i)
                        case 1, ab.x = s.disp.xcen - 10; ab.y = s.disp.ycen + 10; ang = -40; 
                        case 2, ab.x = s.disp.xcen - 10; ab.y = s.disp.ycen + 10; ang = -20; 
                        case 3, ab.x = s.disp.xcen - 10; ab.y = s.disp.ycen - 10; ang = +20; 
                        case 4, ab.x = s.disp.xcen - 10; ab.y = s.disp.ycen - 10; ang = +40; 
                        case 5, ab.x = s.disp.xcen + 10; ab.y = s.disp.ycen - 10; ang = +40; 
                        case 6, ab.x = s.disp.xcen + 10; ab.y = s.disp.ycen - 10; ang = +20; 
                        case 7, ab.x = s.disp.xcen + 10; ab.y = s.disp.ycen + 10; ang = -20; 
                        case 8, ab.x = s.disp.xcen + 10; ab.y = s.disp.ycen + 10; ang = -40; 
                    end
                    % redraw fixation cross and confidence butterfly
                    Screen(w, 'DrawTexture', fix.texture, [], fix.rect);
                    dfi_conf_butterfly_draw(s);
                    % highlight response button
                    ab.resp_angle = ang;
                    dfi_conf_butterfly_fill(s, ab);
                    Screen('DrawTexture', s.disp.win, fix.texture, [], fix.rect);
                    Screen('Flip', w)
                end
            end
            KbEventFlush(s.kb.subj)
        end
    end % end response collection

    % mark end of trial in eyelink log
    if s.el.online
        Eyelink('Message', 'END TRIAL');
        [saccadeOccured] = dfi_checkForSaccades(s);
        if saccadeOccured
            saccCount = saccCount + 1;
            % Feedback about central fixation after a trial
            if strcmp(s.eyelink.fb, 'trial')
                Screen(s.disp.win, 'Flip'); 
                DrawFormattedText(s.disp.win, 'Fixate on the center', 'center', 'center', s.disp.white);
                Screen(s.disp.win, 'Flip'); 
                WaitSecs(1.5); 
                Screen(s.disp.win, 'DrawTexture', s.stim.fix.texture, [], s.stim.fix.rect);
                Screen(s.disp.win, 'Flip'); 
            end
        end
    end
    
    
    %.......................
    %% Staircase response
    %.......................
    if strcmp(s.stair.incl, 'yes')
        % update chain with response
        resp_acc = 0;
        if data.resp(i) > 0
            if (data.resp(i) >= 5 && ismember(trlid(i), [3,6]) && strcmp(data.rkeycfg(i),'AB')) || ...
               (data.resp(i) <  5 && ismember(trlid(i), [3,6]) && strcmp(data.rkeycfg(i),'BA')) || ...
               (data.resp(i) >= 5 && ismember(trlid(i), 8) && strcmp(data.rkeycfg(i),'BA')) || ...
               (data.resp(i) <  5 && ismember(trlid(i), 8) && strcmp(data.rkeycfg(i),'AB'))
                resp_acc = 1;
            end
        end
        % don't count this trial if there were multiple responses
        if data.badRT(i) == 0 && data.resp(i) > 0
            if ismember(trlid(i),3), s.ud{chainID}.FF  = PAL_AMUD_updateUD(s.ud{chainID}.FF,  resp_acc); end;
            if ismember(trlid(i),6), s.ud{chainID}.Fus = PAL_AMUD_updateUD(s.ud{chainID}.Fus, resp_acc); end;
            if ismember(trlid(i),8), s.ud{chainID}.Fis = PAL_AMUD_updateUD(s.ud{chainID}.Fis, resp_acc); end;
            % Did the desired number of reversals occur within the desired
            % number of trials, if specified? We might also burn the first
            % few reversals
            if isfield(s.stair, 'nTrials')
                ntrl = s.stair.nTrials;
                nrvs = s.stair.nRevers;
                burn = s.stair.nBurnin;
                if length(s.ud{chainID}.FF.reversal) > ntrls
                    if sum(s.ud{chainID}.FF.reversal(end-ntrl:end)>0)  < nrvs, s.ud{chainID}.FF.stop  = 0; end 
                    if sum(s.ud{chainID}.FF.reversal>0)  < nrvs + burn, s.ud{chainID}.FF.stop  = 0;        end 
                end
                if length(s.ud{chainID}.Fus.reversal) > ntrls
                    if sum(s.ud{chainID}.Fus.reversal(end-ntrl:end)>0) < nrvs, s.ud{chainID}.Fus.stop = 0; end
                    if sum(s.ud{chainID}.Fus.reversal>0) < nrvs + burn, s.ud{chainID}.Fus.stop = 0;        end
                end
                if length(s.ud{chainID}.Fis.reversal) > ntrls
                    if sum(s.ud{chainID}.Fis.reversal(end-ntrl:end)>0) < nrvs, s.ud{chainID}.Fis.stop = 0; end
                    if sum(s.ud{chainID}.Fis.reversal>0) < nrvs + burn, s.ud{chainID}.Fis.stop = 0;        end
                end
            end
            % when all cues have converged, stop stimulus presentation
            stop = 1;
            for ii = 1:s.stair.nChains
                if s.ud{ii}.FF.stop == 0,  stop = 0; end;
                if s.ud{ii}.Fus.stop == 0, stop = 0; end;
                if s.ud{ii}.Fis.stop == 0, stop = 0; end;
            end
            if stop
                s.stair.fullstop = 1;
                return;
            end;
        end
    end
    
end % end trial loop 

% Feedback about central fixation after a block
if saccCount > 1 && strcmp(s.eyelink.fb, 'block')
    Screen(s.disp.win, 'Flip'); 
    DrawFormattedText(s.disp.win, 'Please try to fixate on the cross in the middle of the screen.', 'center', 'center', s.disp.white);
    Screen(s.disp.win, 'Flip'); 
    WaitSecs(3); 
end


%++++++++++++++++++++
%% 3.)   Finish
%++++++++++++++++++++

% send trigger for end of block
dfi_sendTrigger(s, s.trig.eBlock)

% Empty keyboard queue
KbQueueFlush(s.kb.subj);

% Flip duration:
missm = missm > 0;
fldur = ftsm - vbltm;

if s.fb.flipdur    
    
    disp('Estimated missed flips: ...')
    disp(dataset(missm(:,1), missm(:,2), missm(:,3), missm(:,4), missm(:,5), missm(:,6), ...
                            'VarNames', {'Fix1', 'Fix2', 'V1', 'V1off', 'V2', 'V2off'}))
    disp('Estimated flip duration: ...')
    disp(dataset(fldur(:,1), fldur(:,2), fldur(:,3), fldur(:,4), fldur(:,5), fldur(:,6), ...
                            'VarNames', {'Fix1', 'Fix2', 'V1', 'V1off', 'V2', 'V2off'}))
end

% collect timestamps and responses
time = struct('vbltm',vbltm, 'sotm',sotm, 'fldur', fldur, 'missm', missm);

% prepare trial labels for plotting
lbl= dataset(trlid(1:ntrls), 'VarNames', {'trltypeVA'});

% collect inter event intervals 
idx    = dataset((1:ntrls)', 'VarNames', {'Idx'});
dslack = dataset(repmat(slack,ntrls,1),vbltm(:,2) - vbltm(:,1),  'VarNames', {'ExpSlack', 'RealSlack'});
dV1    = dataset(ITI,                    vbltm(:,3) - vbltm(:,2),  'VarNames', {'ExpV1on', 'V1on'});
dv1dur = dataset(repmat(vdur,ntrls,1), vbltm(:,4) - vbltm(:,3),  'VarNames', {'ExpV1off', 'V1off'});
dA1    = dataset(ITI+ifi/2,             sAudio(:,1) - vbltm(:,2),  'VarNames', {'ExpA1on', 'A1on'});
dsoaA  = dataset(soa,                   sAudio(:,2) - sAudio(:,1), 'VarNames', {'Expsoa', 'Asoa'});
dV2    = dataset(                        vbltm(:,5) - vbltm(:,3),  'VarNames', {'Vsoa'});
dv2dur = dataset(repmat(vdur,ntrls,1), vbltm(:,6) - vbltm(:,5),  'VarNames', {'ExpV2off', 'V2off'});

if s.fb.flipdiagbig
    % print overview of inter-event-intervals
    disp('times between vbl retraces:')
    disp('Slack:')
    disp([idx, lbl, dslack, dataset(dslack.RealSlack-dslack.ExpSlack, 'VarNames', {'Diff'})])
    disp('V1:')
    disp([idx, lbl, dV1, dataset(dV1.V1on-dV1.ExpV1on, 'VarNames', {'Diff'})])
    disp('V1.dur:')
    disp([idx, lbl, dv1dur, dataset(dv1dur.V1off-dv1dur.ExpV1off, 'VarNames', {'Diff'})])
    disp('A1:')
    disp([idx, lbl, dA1, dataset(dA1.A1on-dA1.ExpA1on, 'VarNames', {'Diff'})])
    disp('soa (A):')
    disp([idx, lbl, dsoaA, dataset(dsoaA.Asoa-dsoaA.Expsoa,  'VarNames', {'Diff'})])
    disp('V2:')
    disp([idx, lbl,dsoaA(:,1),dV2, dataset(dV2.Vsoa-dsoaA.Expsoa, 'VarNames', {'Diff'})])
    disp('V2.dur:')
    disp([idx, lbl, dv2dur, dataset(dv2dur.V2off-dv2dur.ExpV2off, 'VarNames', {'Diff'})])
end

% save dataset as backup and to be able to quickly check stimulus timing precision
savepath = fullfile(s.mypath.data,'backup_printout' );
fname = sprintf('printout_subj%s_sess%g_run%g_%s.mat', s.subj.id, s.subj.session, s.runid, datestr(now, 'dd-mmm-yyyy'));
if ~exist(savepath, 'dir'), mkdir(savepath); end;  
if exist(fullfile(savepath,fname), 'file')
    load(fullfile(savepath,fname))
end
printout{b} = [lbl, dslack, dataset(dslack.RealSlack-dslack.ExpSlack, 'VarNames', {'Diff_Slack'}), ...
                    dV1,    dataset(dV1.V1on-dV1.ExpV1on,             'VarNames', {'Diff_V1on'}),  ...
                    dv1dur, dataset(dv1dur.V1off-dv1dur.ExpV1off,     'VarNames', {'Diff_V1off'}), ...
                    dA1,    dataset(dA1.A1on-dA1.ExpA1on,             'VarNames', {'Diff_A1on'}),  ... 
                    dsoaA,  dataset(dsoaA.Asoa-dsoaA.Expsoa,          'VarNames', {'Diff_Asoa'}),  ...
                    dV2,    dataset(dV2.Vsoa-dsoaA.Expsoa,            'VarNames', {'Diff_V2on'}),  ...
                    dv2dur, dataset(dv2dur.V2off-dv2dur.ExpV2off,     'VarNames', {'Diff_V2off'}), ...
               ];
p = printout{b};
printout_small{b} = dataset(p.trltypeVA, p.Diff_Slack, p.Diff_V1on, p.Diff_V1off, ...
                     p.Diff_A1on, p.Diff_Asoa, p.Diff_V2on, p.Diff_V2off, ...
                     'VarNames', {'trlid', 'Diff_Slack', 'Diff_V1on', 'Diff_V1off' ...
                          'Diff_A1on', 'Diff_Asoa', 'Diff_V2on', 'Diff_V2off'});
save(fullfile(savepath, fname), 'printout', 'printout_small')

if s.fb.flipdiagsmall
    disp('small printout:')
    disp(printout_small{b})
end


%+++++++++++++++++++++++++++++
%% 4.)   Nested functions
%+++++++++++++++++++++++++++++

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

% 
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






% 
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    