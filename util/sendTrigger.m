function sendTrigger( s, trig )
%sendTrigger 
% uses a LabJack device to send triggers to an external hardware, like an
% EEG setup. Needs labjack object (lj) and the trigger value to be send (trig). 

    % Send trigger now. 
    sendNow = 1;
    prepareStrobe(s.lj, trig, [], sendNow);
    
end
% eof

