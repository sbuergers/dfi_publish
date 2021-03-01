function [time, data] = dfi_runSegment_test(s, data, b)
%[s, d] = dfi_runSegment(s, data)
% starts the sequence of stimuli of one block of the VE paradigm (see
% dfi_run_experiment for more detailed information on the paradigm)
%
% DFI:
%                     A            A
% |-------------------V-----------(V)----------------R---------------|-------
%    ITI(650-1950ms)       soa        reswin(2000ms)    slack(800ms)    ITI
%
% FE:
%                                 
% |-------------------V------------V----------------R---------------|-------
%    ITI(650-1950ms)       soa       reswin(2000ms)    slack(800ms)    ITI
%
% -----
% adapted from Agoston Mihalic, 
% last updated, April 2015
%

% -----------------------------------------
%              some setting up
% -----------------------------------------

% initialize variables
w               = s.disp.win;
fix             = s.stim.fix;
sAudio          = zeros(s.ntrls,2);
trlid           = data.trlid;
trgid           = trlid;

% generate ITIs (mean +- std)
plusminus = Shuffle([ones(s.ntrls/2,1); -ones(s.ntrls/2,1)]);
ITI       = s.ITI(1) + plusminus .* rand(s.ntrls, 1) * s.ITI(2);

% segment soas
soa = data.soa(1:s.ntrls);

% setup flip-times matrix (time x nflips)
% VBLTimestamp StimulusOnsetTime FlipTimestamp Missed (Beampos)
vbltm = zeros(s.ntrls, 6); 
sotm  = zeros(s.ntrls, 6);
ftsm  = zeros(s.ntrls, 6);
missm = zeros(s.ntrls, 6);

% for all intervals get exact multiples of ifi - cutting off the remainder,
% but keeping at least one ifi, as 0 is nonsensical for future times.
ifi   = s.disp.ifi;
slack =                      s.slack - mod(s.slack, ifi);
ITI   =                          ITI - mod(ITI, ifi);
resw  = repmat(s.reswin, s.ntrls, 1) - mod(s.reswin, ifi);
soa   =                          soa - mod(soa, ifi);
soa(soa == 0) = ifi;
vdur  = s.stim.vdur - mod(s.stim.vdur, ifi);
if vdur == 0, vdur = ifi; end;  

% ----------------------------------------
%          Stimulus presentation
% ----------------------------------------

% buffer sound
PsychPortAudio('FillBuffer', s.haudio, s.stim.aud.loc);

% play a silent beep to avoid timing issues with the first sound stimulus
oldVolume = PsychPortAudio('Volume', s.haudio, 0);
PsychPortAudio('Start', s.haudio, 1);  
PsychPortAudio('Stop', s.haudio, 1);
PsychPortAudio('Volume', s.haudio, oldVolume);

% give the buffer some time
WaitSecs(1)

% trial loop
for i=1:s.ntrls
    
    % Abort?
    [keydown,~,keycode] = KbCheck(s.kb.op);
    if keydown && keycode(s.key.quit)
        return
    end
    
    % Fixation cross (idle)
    if i==1
        [vbltm(i,1), sotm(i,1), ftsm(i,1), missm(i,1)] = Screen(w, 'Flip');
    else
        [vbltm(i,1), sotm(i,1), ftsm(i,1), missm(i,1)] = Screen(w, 'Flip', stimStopTime+resw(i-1)-ifi/2);
    end
    
    % Buffer V1
    Screen(w, 'DrawTexture', s.stim.vis.texture, [], s.stim.vis.rect);                            % V1buff
    
    % Display V1/A1
%     sAudio(i,1) = PsychPortAudio('Start', s.haudio, 1, vbltm(i,1)+ITI(1)+ifi/2, waitForSpeakers); % A1 DFI 2 and 3 and 5
    sAudio(i,1) = PsychPortAudio('Start', s.haudio, 1, vbltm(i,1)+ITI(i)+ifi/2+s.stim.vposdelay); % A1 DFI 6

    [vbltm(i,3), sotm(i,3), ftsm(i,3), missm(i,3)] = Screen(w, 'Flip', vbltm(i,1)+ITI(i)-ifi/2);  % V1
%     sAudio(i,1) = PsychPortAudio('Start', s.haudio, 1, vbltm(i,3)+ifi/2, waitForSpeakers);        % A1 DFI 1
%     sAudio(i,1) = PsychPortAudio('Start', s.haudio, 1, vbltm(i,3)+ifi/1.5, waitForSpeakers);        % A1 DFI 4
    dfi_sendTrigger(s, s.trig.sAV + trgid(i))                                                      % Tr.AV
    [vbltm(i,4), sotm(i,4), ftsm(i,4), missm(i,4)] = Screen(w, 'Flip', vbltm(i,3)+vdur-ifi/2);    % V1off
   
    stimStopTime = vbltm(i,4);
    
    % Response
    while data.resp(i)==0 && GetSecs < stimStopTime + resw(i)-ifi
        [keydown, tclick, keycode] = KbCheck(s.kb.subj);
        keyid = find(keycode);
        if keydown && length(keyid) == 1 % multiple responses are excluded
            if keycode(s.key.quit); 
            elseif ismember(keyid, s.key.code)
                data.RT(i)   = tclick - stimStopTime;
                data.resp(i) = find(keyid==s.key.code);
                if data.resp(i) == 1
                    dfi_sendTrigger(s, s.trig.rOne + trgid(i))
                else
                    dfi_sendTrigger(s, s.trig.rTwo + trgid(i))
                end
                resw(i) = data.RT(i) - mod(data.RT(i), ifi);
            end
        end 
    end
end % end trial loop 

% ---------------------------------------
%                 finish
% ---------------------------------------

% Flip diagnostics:
missm = missm > 0;
disp('Estimated missed flips: ...')
disp(dataset(missm(:,1), missm(:,2), missm(:,3), missm(:,4), missm(:,5), missm(:,6), ...
                        'VarNames', {'Fix1', 'Fix2', 'V1', 'V1off', 'V2', 'V2off'}))

disp('Estimated flip duration: ...')
fldur = ftsm - vbltm;
disp(dataset(fldur(:,1), fldur(:,2), fldur(:,3), fldur(:,4), fldur(:,5), fldur(:,6), ...
                        'VarNames', {'Fix1', 'Fix2', 'V1', 'V1off', 'V2', 'V2off'}))
                    
% collect timestamps and responses
time = struct('vbltm',vbltm, 'sotm',sotm, 'fldur', fldur, 'missm', missm);

% prepare trial labels for plotting
lbl= dataset(trlid(1:s.ntrls), 'VarNames', {'trltypeVA'});

% print overview of inter-event-intervals
idx = dataset((1:s.ntrls)', 'VarNames', {'Idx'});
disp('V1:')
dV1 = dataset(ITI, vbltm(:,3) - vbltm(:,1), 'VarNames', {'ExpV1on', 'V1on'});
disp([idx, lbl, dV1, dataset(dV1.V1on-dV1.ExpV1on, 'VarNames', {'Diff'})])
disp('V1.dur:')
dv1dur = dataset(repmat(vdur,s.ntrls,1), vbltm(:,4) - vbltm(:,3), 'VarNames', {'ExpV1off', 'V1off'});
disp([idx, lbl, dv1dur, dataset(dv1dur.V1off-dv1dur.ExpV1off, 'VarNames', {'Diff'})])
disp('A1:')
dA1 = dataset(ITI+ifi/2, sAudio(:,1) - vbltm(:,1), 'VarNames', {'ExpA1on', 'A1on'});
disp([idx, lbl, dA1, dataset(dA1.A1on-dA1.ExpA1on, 'VarNames', {'Diff'})])



function dfi_sendTrigger( s, trig )
%sendTrigger 
% uses a LabJack device to send triggers to an external hardware, like an
% EEG setup. Needs labjack object (lj) and the trigger value to be send (trig). 

    % Send trigger now. 
    if s.ljPresent
%         sendNow = 1;
%         prepareStrobe(s.lj, trig, [], sendNow);
        prepareStrobe(s.lj, trig);
        strobeWord(s.lj);
    end
end
end % eof
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    