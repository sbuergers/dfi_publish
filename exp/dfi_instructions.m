function dfi_instructions(task, s)
% diplays instructions for dfi paradigm
%
% task =
% 'expStart'     : Instructions about the experiment as a whole
% 'expStart2'    : Continued initial instructions
% 'startPractice': Practice coming up
% 'endPractice'  : Practice is over
% 'numFlashes'   : pre-block, attend to flashes
% 'numTones'     : pre-block, attend to tones
% 'endBlock'     : post-block feedback
% 'wrongKeys'    : Last block you pressed the wrong keys
% 'forcedBreak'  : Upcoming is a forced break of ... minutes
% 'endLastBlock' : Farewell 
%
% ----------------------------------------
% adapted from Agoston Mihalik, 
% last updated, June 2016
%
    Screen('TextSize', s.disp.win, 16)

    WaitSecs(0.5); % a bit of extra time
    Screen('FillRect', s.disp.win, s.disp.grey);
    
    % Keyboard response:
    kb_resp_instr_2IFC = sprintf('Use the key "%s" for the 1st interval and "%s" for the 2nd interval.', s.key.resp{:});
    kb_resp_instr_2IFC_AB = sprintf('Use the key "%s" for the 1st interval and "%s" for the 2nd interval.', s.key.resp{:});
    kb_resp_instr_2IFC_BA = sprintf('Use the key "%s" for the 2nd interval and "%s" for the 1st interval.', s.key.resp{:});
    % Mouse response:
    ms_resp_instr_2IFC = sprintf('Click LEFT for the 1st and RIGHT for the 2nd interval.');
    
    % ---------------------------------------------------------
    %                           yesno
    % ---------------------------------------------------------
    if strcmp(s.paradigm, 'yesno') || strcmp(s.paradigm, 'YN_threshold')  
        switch task
            case 'expStart' % Experiment explanation at the very start
                line{1} = 'Welcome. In the following experiment you will be presented with one or two flashes,';
                line{2} = 'sometimes accompanied by one or two tones. Your task is simply to report how many flashes you perceived.';
                line{3} = 'Do not try to reason about how many flashes were presented, we are interested in you perception!';
                line{4} = 'you are reporting one percept more often than another, that is perfectly fine! If you';
                line{5} = 'perceive no flash then do not respond. Both sounds and flashes are always presented on the left side.';
                line{6} = '';
                line{7} = 'During the task, please keep your head resting on the chin-rest, and your eyes fixated on the cross';
                line{8} = 'in the middle of the screen. There will be several breaks during the experiment. Feel free to use that';
                line{9} = 'time to move around a bit, loosen your neck and shoulders and relax your eyes.';
                line{10} = '';
                line{11} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n\n%s\n\n%s\n\n%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
            case 'startPractice' % practice starts now
                line{1} = sprintf('You will now be presented with %g practice trials to make you accustomed to the task.', s.pract.ntrls);
                line{2} = sprintf('');
                line{3} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
            case 'noBlinks' % don't blink in critical period    
                line{1} = 'Please, try not to blink shortly before, during and after stimulus presentation.'; 
                line{2} = 'A good strategy is to blink after the response on each trial.';
                line{3} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
            case 'numFlashes' % before a block      
                if strcmp(s.paradigm, 'yesno')
                    line{1} = 'How many flashes do you see?'; 
                    if strcmp(s.keyswitch, 'AB')
                        line{2} = sprintf('Use the key "%s" for 1 flash and "%s" for 2 flashes.', s.key.resp{:});
                    else
                        line{2} = sprintf('Use the key "%s" for 2 flashes and "%s" for 1 flash.', s.key.resp{:});
                    end
                    line{3} = 'Press any button to continue!';
                    DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
                elseif strcmp(s.paradigm, 'YN_threshold')
                    line{1} = 'How many flashes do you see?'; 
                    if strcmp(s.keyswitch, 'AB')
                        line{2} = sprintf('Use the keys "%s", "%s", "%s", and "%s" for 1 flash and "%s", "%s", "%s", and "%s" for 2 flashes.', s.key.resp{:});
                    else
                        line{2} = sprintf('Use the keys "%s", "%s", "%s", and "%s" for 2 flashes and "%s", "%s", "%s", and "%s" for 1 flash.', s.key.resp{:});
                    end    
                    line{3} = 'Press any button to continue!';
                    DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
                end
            case 'numTones'   % before a block
                line{1} = 'How many sounds do you hear?'; 
                line{2} = sprintf('Use the keys "%s" for 1 and "%s" for 2 tones.', s.key.resp{:});
                line{3} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
            case 'expStart2' % Continued experiment explanation 
                if ~strcmp(s.cond, 'multiAttAud')
                    line{1} = 'Note that in some cases you might miss the flash(es), because you blinked at the presentation time.';
                    line{2} = 'If that happens, simply wait for the next trial and do not respond.';
                    line{3} = 'There is a maximum response window of 1.5 seconds. Occasionally there will be trials with no stimuli.';
                    line{4} = 'Press any button to continue!';
                    if strcmp(s.paradigm, 'YN_threshold')
                        line{3} = 'There is a maximum response window of 2.5 seconds.';
                    end
                else
                    line{1} = 'Note that this time you should attend to the sounds and ignore the flashes!';
                    line{2} = 'Nevertheless, fixate on the cross and do NOT close your eyes.';
                    line{3} = 'There is a maximum response window of 1.5 seconds.';
                    line{4} = 'Press any button to continue!';
                end
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
        end % switch yesno
    % ---------------------------------------------------------
    %                           2AFC
    % ---------------------------------------------------------
    elseif strcmp(s.paradigm, '2AFC')
        switch task
            case 'expStart' % Experiment explanation at the very start
                line{1} = 'Welcome. In the following experiment you will be shown two intervals of one or two flashes,';
                line{2} = 'sometimes accompanied by concurrent tones. You should determine which interval contained two flashes.';
                line{3} = 'Both sounds and flashes are always presented on the left side.';
                line{4} = 'During the task, please keep your head resting on the chin-rest, and your eyes fixated on the cross';
                line{5} = 'in the middle of the screen. There will be several breaks during the experiment. Feel free to use that';
                line{6} = 'time to move around a bit, loosen your neck and shoulders and relax your eyes.';
                line{7} = '';
                line{8} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n%s\n\n\n%s\n\n%s\n\n%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
            case 'numFlashes' % before a block
                line{1} = 'Which interval contained 2 flashes?';
                if strcmp(s.respdev, 'keyboard')
                    if strcmp(s.keyswitch, 'AB')
                        line{2} = kb_resp_instr_2IFC_AB;
                    elseif strcmp(s.keyswitch, 'BA')
                        line{2} = kb_resp_instr_2IFC_BA;
                    end
                elseif strcmp(s.respdev, 'mouse')
                    line{2} = ms_resp_instr_2IFC;
                end
                line{3} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
            case 'exampleTrials' 
                line{1} = 'You will now be shown a few examples of either'; 
                line{2} = '<flash> <flash> .... <flash>; or';
                line{3} = '<flash> .... <flash> <flash> trials. Try to find the correct interval';
                line{4} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
            case 'exampleFeedbackCorr' 
                line{1} = 'Correct.'; 
                line{2} = '';
                line{3} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
            case 'exampleFeedbackIncorr' 
                line{1} = 'Incorrect.'; 
                line{2} = '';
                line{3} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
            case 'startPractice' % practice starts now
                line{1} = sprintf('You will now be presented with %g practice trials. There will be no feedback anymore', s.pract.ntrls);
                line{2} = sprintf('and trials are no longer self paced so it might feel a bit fast at first.');
                line{3} = sprintf('You will also notice that the 2 flash interval is much harder to see in many cases.');
                line{4} = sprintf('Simply make your best guess when unsure. To make it easier to notice when a new trial has started');
                line{5} = sprintf('the fixation cross rotates at trial onset. Finally, when you are sure there were two flashes in the');
                line{6} = sprintf('first interval, please WAIT with your response until you have seen the other interval as well.');
                line{7} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
            case 'numTones'   % before a block
                line{1} = 'Which interval contained 2 tones?'; 
                if strcmp(s.respdev, 'keyboard')
                    line{2} = kb_resp_instr_2IFC;
                elseif strcmp(s.respdev, 'mouse')
                    line{2} = ms_resp_instr_2IFC;
                end
                line{3} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
        end % switch 2AFC
    end % end if yesno
    
    % ---------------------------------------------------------  
    %                   Condition unspecific
    % ---------------------------------------------------------
    switch task
        case 'confidenceInstruction_prePractice'
            if strcmp(s.respdev, 'keyboard')
                line{1} = 'For this task, in addition to your percept you indicate your confidence in your percept.';
                line{2} = 'Your confidence is represented by your response fingers, with increasing confidence from pinkie to index finger.';
                line{3} = 'For the upcoming practice you will see your responses mapped on the screen to make this clearer. Later no feedback is given.';
                line{4} = 'Press any button to continue!';
                DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
            elseif strcmp(s.respdev, 'mouse')
                % to do...
            end
        case 'response_finger_position'
            line{1} = 'To respond, please put your hands in the 10-finger starting position.';
            line{2} = 'That is asdf (left hand), jkl; (right hand).';
            line{3} = 'Before each block you will be instructed, which keys correspond to which answer.';
            line{4} = 'Press any button to continue!';
            DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
        case 'endPractice'   % practice is over
            line{1} = 'This is the end of the practice period.';
            line{2} = 'Now the experimental part is going to start.';
            line{3} = 'If you have any remaining questions, please ask them now.';
            line{4} = 'Press any button to continue!';
            DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
        case 'endPracticeBlock1'   % practice is over
            line{1} = 'This is the end of the first out of two practice blocks.';
            line{2} = 'Please, call the experimenter before continuing.';
            line{3} = '';
            line{4} = 'Press any button to continue!';
            DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n%s\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
        case 'endBlock'      % after block feedback
            line{1} = sprintf('This is the end of block %g out of %g. Feel free to take a short break.', [s.curbl, s.blocks.N]);
            line{3} = 'Press any button to continue!';
            DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
        case 'wrongKeys'     
            line{1} = sprintf('Unfortunately, it seems that you have been pressing the incorrect response keys.');
            line{3} = 'Press any button to continue!';
            DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
        case 'forcedBreak'
            line{1} = sprintf('This is the end of block %g out of %g, and a forced break.', [s.curbl, s.blocks.N]);
            line{2} = sprintf('You will not be able to continue for the next %g seconds.', s.breakdur);
            line{3} = 'Press any button to continue!';
            DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
        case 'calibrationComplete'
            line{1} = sprintf('Calibration complete.');
            line{3} = 'Press any button to continue!';
            DrawFormattedText(s.disp.win, sprintf('%s\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
        case 'driftCheck'
            line{1} = sprintf('Drift check: Please focus on the fixation cross.');
            line{3} = 'Press any button to continue!';
            DrawFormattedText(s.disp.win, sprintf('%s\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
        case 'endLastBlock'   % Experiment is over (for now)!
            line{1} = 'This is the end of the current run.';
            line{2} = 'Well done!';
            line{3} = 'Press any button to continue!';
            DrawFormattedText(s.disp.win, sprintf('%s\n\n%s\n\n\n\n%s', line{:}), 'center', 'center', s.disp.white);
    end % switch 
    
    % display on screen
    Screen(s.disp.win, 'Flip');
    if strcmp(task, 'forcedBreak')
        WaitSecs(s.breakdur); 
    end;
    KbWait(s.kb.subj); % wait for response  
    if s.disp.incolor
        Screen(s.disp.win, 'FillRect', s.disp.green);
    else
        Screen(s.disp.win, 'FillRect', s.disp.grey);
    end
    WaitSecs(0.5);
    Screen(s.disp.win, 'Flip');

end % end dfi_instructions
% eof











