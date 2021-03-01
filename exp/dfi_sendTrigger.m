function dfi_sendTrigger( s, trig )
%sendTrigger 
% uses a LabJack device to send triggers to an external hardware, like an
% EEG setup. Needs labjack object (s.lj) and the trigger value to be send (trig). 
    % Send trigger now. 
    if s.ljPresent
        prepareStrobe(s.lj, trig);
        strobeWord(s.lj);
    end
    % Send eyelink trigger (pretty accurate timekeeping, 1-2 ms delay)
    if s.el.online
        Eyelink('Message', num2str(trig));
    end
end