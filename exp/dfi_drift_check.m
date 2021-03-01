function driftDistance = dfi_drift_check(s, tfix, ttimeout)
% estimates the discrepancy between central fixation and actual gaze, t is
% the time taken to estimate fixation coordinates. If a saccade occurs
% within that time the function restarts the estimation process. Time out
% is by default after 10 seconds. 


% default values for input parameters
if ~exist('tfix', 'var')
    tfix = 1;
end

if ~exist('ttimeout', 'var')
    ttimeout = 10;
end    
    
    % make sure central fixation is maintained for at least a second and
    % there are no major saccades during that time
    centralFixation = 0;
    totalTime = 0;
    pos = [];
    tstart = GetSecs;
    trun   = GetSecs;
    
    while ~centralFixation
        
        % initialize some variables
        drained = 0;
        gazeavg = [];
        saccend = [];
        
        % only works if we have an eye to track so
        try
            while ~drained % get all queued data

                [~, events, drained] = Eyelink('GetQueuedData', s.el.LEFT_EYE);

                if ~isempty(events) && any(ismember(events(2,:), [6,9]))

                    %% do we have some undesirable saccades?
                    % get gaze positions (fixation updates and saccades)
                    if sum((events(2,:) == 9)) > 0
                        gazeavg(:,1) = events(19, (events(2,:) == 9));
                        gazeavg(:,2) = events(20, (events(2,:) == 9));                    
                    end

                    if sum((events(2,:) == 6)) > 0
                        saccend(:,1) = events(14, (events(2,:) == 6));
                        saccend(:,2) = events(15, (events(2,:) == 6));                    
                    end

                    pos = [pos; gazeavg; saccend];

                    allowedDev = dfi_deg2pix([s.eyelink.driftdeg+1 0], s.disp);

                    % maximally allowed saccade distance (circle around fix)
                    if any(  sqrt((pos(:,1) - s.disp.xcen).^2 + (pos(:,2) - s.disp.ycen).^2) > allowedDev(1)  ) 
                        pos = [];
                    end
                end

            end % while drain queue
        end
        
        % continue doing this until fixation is maintained for time tfix
        elapsedTime = GetSecs - trun;
        
        if ~isempty(pos) 
            centralFixation = 1; 
        else
            trun = GetSecs;
        end
                
        if elapsedTime < tfix
            centralFixation = 0; 
        end;
        
        % abort if too much time has elapsed       
        if GetSecs - tstart > ttimeout
            driftDistance = 999; 
            return;
        end
        
    end % while central fixation
    
    % calculate distance
    driftDistance = mean( sqrt( (pos(:,1) - s.disp.xcen).^2 + (pos(:,2) - s.disp.ycen).^2 ) );
    fprintf('DriftDistance is %d\n\n', driftDistance);
    
    
return % end fun

% eof














