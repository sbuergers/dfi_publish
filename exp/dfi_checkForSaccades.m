function [saccadeOccured, events] = dfi_checkForSaccades(s)
% #define STARTBLINK 3    // pupil disappeared, time only
% #define ENDBLINK   4    // pupil reappeared, duration data
% #define STARTSACC  5    // start of saccade, time only
% #define ENDSACC    6    // end of saccade, summary data
% #define STARTFIX   7    // start of fixation, time only
% #define ENDFIX     8    // end of fixation, summary data
% #define FIXUPDATE  9    // update within fixation, summary data for interval
% 
% #define MESSAGEEVENT 24  // user-definable text: IMESSAGE structure
% 
% #define BUTTONEVENT  25  // button state change:  IOEVENT structure
% #define INPUTEVENT   28  // change of input port: IOEVENT structure

    drained = 0;
    saccadeOccured = 0;
    numDrains = 0;
    gazeavg = [];
    saccend = [];
    while ~drained
        [~, events, drained] = Eyelink('GetQueuedData', s.el.LEFT_EYE);
        
        if ~isempty(events)
            if any(ismember(events(2,:), [6,9])) % end sacc, fixupdate (all recorded gaze positions)
                
                % get gaze positions (fixation updates and saccades)
                if sum((events(2,:) == 9)) > 0
                    gazeavg(:,1) = events(19, (events(2,:) == 9));
                    gazeavg(:,2) = events(20, (events(2,:) == 9));                    
                end
                    
                if sum((events(2,:) == 6)) > 0
                    saccend(:,1) = events(14, (events(2,:) == 6));
                    saccend(:,2) = events(15, (events(2,:) == 6));                    
                end
                
                pos = [gazeavg; saccend];

                allowedDev = dfi_deg2pix([s.eyelink.fixdeg 0], s.disp);

%                 % maximally allowed saccade distance (circle around fix)
%                 if any(  sqrt((pos(:,1) - s.disp.xcen).^2 + (pos(:,2) - s.disp.ycen).^2) > allowedDev(1)  ) 
%                     saccadeOccured = 1;     
%                     fprintf('\nSaccade deteced...!\n');
%                     pos
%                 end
                
                % maximally allowed saccade distance (focus on x)
                idx = find(  sqrt((pos(:,1) - s.disp.xcen).^2 + (pos(:,2) - s.disp.ycen).^2) > allowedDev(1)  );
                if any(  (pos(idx,1) - s.disp.xcen).^2 - (pos(idx,2) - s.disp.ycen).^2 > 0  )
                    saccadeOccured = 1;     
                    fprintf('\nHorizontal saccade detected...!\n');
                    pos
                end
                
            end
        end
        
        numDrains = numDrains + 1;
%         fprintf('Number of drains: %d\n\n', numDrains);
    end
end
















