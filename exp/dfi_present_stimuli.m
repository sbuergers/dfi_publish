function [s, data] = dfi_present_stimuli(s, debug)
%[s, data] = dfi_present_stimuli(s, debug)
% takes parameters from dfi_run_experiment and executes an audiovisual
% integration task. It contains the block-loop. Calls the functions
%
% dfi_setupPTB
% dfi_prepFixCross
% dfi_prepViStim
% dfi_prepAuStim
% dfi_setupKeyboard
% dfi_instructions
% dfi_runSegment
%
% ------
% adapted from Agoston Mihalik, 
% last updated, April 2015
%           

% without sound only condition (see bottom for with sound only condition)
if s.practice
    % pre-practice (sort of part of the instructions) for 2IFC task
    if strcmp(s.paradigm, '2AFC')
        % Introduce the task and make a few practice trials with feedback:
        dfi_instructions('exampleTrials', s); 
        if ~strcmp(s.cond, 'multiAttAud')
            dfi_instructions('numFlashes', s);
        else 
            dfi_instruction('numTones', s);
        end
        datatemp = dfi_generate_data(s, 1); 
        ppdata   = datatemp(1:6,:);
        ppdata.soa = [0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25];
        ppdata.trlid = Shuffle([2,3,5,6,8,9,2,3,5,6,8,9]');
        spract   = s;
        spract.practicefeedback = 1;
        spract.ntrls = 12; 
        spract.stim     = dfi_prepViStim(s, ppdata);
        spract.stim     = dfi_prepAuStim(spract.stim);
        spract.stim.fix = dfi_prepFixCross(spract);
        [data, abort]   = dfi_runSegment_2AFC(spract, ppdata, 1);
    end
        
    % generate practice data
    sp = s;
    sp.blocks.N   = s.pract.blocks.N; 
    sp.blocks.Nv  = s.pract.blocks.Nv;
    sp.blocks.Nav = s.pract.blocks.Nav;
    sp.trlid = [2 3 2 3]; pdata{1} = dfi_generate_data(sp, 1); % flashes
    sp.trlid = [5 6 8 9]; pdata{2} = dfi_generate_data(sp, 1); % both
    % randomly show Vis, Aud and multimodal practice blocks
    if strcmp(s.paradigm, 'YN_threshold')
        dfi_instructions('confidenceInstruction_prePractice', s)
    end
    dfi_instructions('startPractice', s) 
    loopNum = 0; 
    for practCond = Shuffle([1 2])
        loopNum = loopNum + 1;
        if ~strcmp(s.cond, 'multiAttAud')
            dfi_instructions('numFlashes', s);
        else 
            dfi_instructions('numTones', s);
        end           
        spract = sp; spract.ntrls = s.pract.ntrls;
        % prepare stimuli
        spract.stim     = dfi_prepViStim(spract, pdata{practCond});
        spract.stim     = dfi_prepAuStim(spract.stim); 
        spract.stim.fix = dfi_prepFixCross(spract);
        if strcmp(s.paradigm, 'yesno') || strcmp(s.paradigm, 'YN_threshold')
            [~, data, abort] = dfi_runSegment(spract, pdata{practCond}, 1);
        elseif strcmp(s.paradigm, '2AFC')
            [data, abort] = dfi_runSegment_2AFC(spract, pdata{practCond}, 1); 
        end
        if loopNum == 1
           dfi_instructions('endPracticeBlock1', s) 
        end
            
        % did we press the abort key q?
        if abort
            return
        end
    end 
    
    % turn wifi on
    tic
    system('networksetup -setairportpower en1 on');

    % Upload practice data
    fprintf('Uploading data to google drive after block .....\n');
    try
        DrawFormattedText(s.disp.win, 'Uploading data to google drive. This can take a few seconds.....', 'center', 'center', s.disp.white);
        Screen('Flip', s.disp.win);
        upload_path = fullfile('/Users/sxb1173/Google Drive/SIFI_interim_data', s.subj.id, s.paradigm);
        if ~exist(upload_path, 'dir'), mkdir(upload_path); end;
        save(fullfile(upload_path, sprintf('data_subj%s_%s_sess%d_run%d%_practice',s.subj.id, s.paradigm, s.subj.session, s.runid)), 'data');
        WaitSecs(10);
        fprintf('.....successful!\n');
    catch
        fprintf('.....unsuccessful!\n');
    end
    toc

    % turn wifi off again
    system('networksetup -setairportpower en1 off')
    
    dfi_instructions('endPractice', s)
    s.practice = 0;
    return
end

% -----------------------------------
%       experimental trials
% -----------------------------------
if ~s.practice
    % generate data
    data = dfi_generate_data(s, 0); 
    
    % STAIRCASE add chainID column in data if staircase is included
    if strcmp(s.stair.incl, 'yes')
        data.chainid = zeros(size(data,1),1);
    end

    time = [];
    for b=1:s.blocks.N

        % keep track of block for instructions
        s.curbl = b;

        % stimulus data of current block
        dbl = data(1+(b-1)*s.ntrls:b*s.ntrls,:);

        % prepare stimuli
        s.stim     = dfi_prepViStim(s, dbl); 
        s.stim     = dfi_prepAuStim(s.stim); 
        s.stim.fix = dfi_prepFixCross(s);
        
        % Eyelink: this is the start of a new block, do we have to recalibrate?
        if s.el.online
            % do custom drift check, for this we need link data
            Eyelink('Message', 'DRIFT CHECK START, BLOCK %d', b); % send 'BLOCKID'
            Eyelink('StartRecording', [], [], [], 1);
            % present fixation stimulus
            dfi_instructions('driftCheck', s)
            Screen(s.disp.win, 'DrawTexture', s.stim.fix.texture, [], s.stim.fix.rect);
            Screen(s.disp.win, 'Flip');
            WaitSecs(0.5)
            % get average fixation in a 1 second window
            allowedDrift = dfi_deg2pix([s.eyelink.driftdeg 0], s.disp);
            drift = dfi_drift_check(s, 0.5, 10);
            % force recalibration after practice
            if b == 1, drift = 999; end;
            if  drift > allowedDrift(1)
                EyelinkDoTrackerSetup(s.el);
                EyelinkDoDriftCorrection(s.el);
                dfi_instructions('calibrationComplete', s)
            end
            Eyelink('Message', 'DRIFT CHECK FINISH, BLOCK %d', b); % send 'BLOCKID'
            Eyelink('StopRecording');
        end

        % instruct
        if any(ismember(s.stim.trlid, [4,7])) || strcmp(s.cond, 'multiAttAud')
            dfi_instructions('numTones', s) 
        else
            dfi_instructions('numFlashes', s) 
        end 
        
        % start eyelink recording
        if s.el.online
            Eyelink('Message', 'BLOCKID %d', b); % send 'BLOCKID'
            Eyelink('StartRecording'); 
            WaitSecs(0.01);
            Eyelink('Message', 'SYNCTIME');

            % are we recording?
            if Eyelink('CheckRecording'); % Returns 0 if recording in progress
                Eyelink('Shutdown');
                s.el.online = 0;
            end
        end
        
        % turn wifi off during task performance
        system('networksetup -setairportpower en1 off')

        % run segments
        if strcmp(s.stim.vis.type, 'test')
            [tblock, dblock] = dfi_runSegment_test(s, dbl, b);
        elseif strcmp(s.paradigm, 'yesno') || strcmp(s.paradigm, 'YN_threshold')
            [tblock, dblock, abort, s] = dfi_runSegment(s, dbl, b);
        elseif strcmp(s.paradigm, '2AFC')
            [dblock, abort] = dfi_runSegment_2AFC(s, dbl, b);
            tblock = [];
        elseif strcmp(s.paradigm, '2AFCsoundTest')
            dblock = dfi_runSegment_2AFC_soundTest(s, dbl, b);
            tblock = [];
        end
        
        % close all off screen windows
        Screen('Close');

        % clear java heap memory
        disp('Clearing Java heap memory...');
        jheapcl;

        % collect data
        data(1+(b-1)*s.ntrls:b*s.ntrls,:) = dblock;
        time                              = [time; tblock];
        
        % did we press the abort key q?
        if abort
            return
        end
        
        % turn wifi on again
        tic
        system('networksetup -setairportpower en1 on');
        
        % Eyelink, stop listening!
        if s.el.online
            if Eyelink('IsConnected')
                WaitSecs(0.015); % capture last events
                Eyelink('StopRecording');
            elseif s.el.online
                warning('Lost EyeLink connection!')
            end
        end
        
        % Upload block data
        fprintf('Uploading data to google drive after block .....\n');
        try
            DrawFormattedText(s.disp.win, 'Uploading data to google drive after block. This can take a few seconds.....', 'center', 'center', s.disp.white); 
            Screen('Flip', s.disp.win); 
            upload_path = fullfile('/Users/sxb1173/Google Drive/SIFI_interim_data', s.subj.id, s.paradigm);
            if ~exist(upload_path, 'dir'), mkdir(upload_path); end;
            save(fullfile(upload_path, sprintf('data_subj%s_%s_run%d_block%d',s.subj.id, s.paradigm, s.runid, b)), 'dblock', 's'); 
            WaitSecs(10)
            fprintf('.....successful!\n');
        catch
            fprintf('.....unsuccessful!\n');
        end
        toc
        
        % STAIRCASE: if all chains converged, shutdown
        if strcmp(s.stair.incl, 'yes') && s.stair.fullstop == 1
            return 
        end;
        
        % Did you press the wrong keys?
        if sum(dblock.RT == 0) > numel(dblock.RT) / 4
            % in the staircase non-responses are expected once chains start
            % to converge
            if ~strcmp(s.stair.incl, 'yes')
                dfi_instructions('wrongKeys', s)
            end
        end

        % show instructions
        if b < s.blocks.N
            if mod(b, s.breakint) > 0
                dfi_instructions('endBlock', s)
            else
                dfi_instructions('forcedBreak', s)
            end
        end
    end % end segment loop

    % add time to setup struct
    s.time = time;
    return
    
end % no practice

% eof






