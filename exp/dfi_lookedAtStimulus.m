function [peakedAtStimulus, events] = dfi_lookedAtStimulus(s)
% gives TRUE for peakedAtStimulus if central fixation wanders off more than
% 's.eyelink.fixdegrad' degrees to the left.

    [~, events] = Eyelink('GetQueuedData', s.el.LEFT_EYE);
    peakedAtStimulus = 0;
    if ~isempty(events)
        if any(ismember(events(2,:), 6))
            % get start and end coordinates
            % in pixels
            spos(:,1) = events(9,  (events(2,:) == 6));
            spos(:,2) = events(10, (events(2,:) == 6));
            epos(:,1) = events(14, (events(2,:) == 6));
            epos(:,2) = events(15, (events(2,:) == 6));
             
            MaxErrLeft  = dfi_deg2pix([s.eyelink.fixdegrad, 0], s.disp);
            
            if any(epos(:,1) - spos(:,1)) > MaxErrLeft
                peakedAtStimulus = 1;
            end
        end
    end
% lets assume that it will always drain completely if I call it for every
% trial (or at least will include all relevant last trial data)
end