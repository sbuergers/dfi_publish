function [i_NoEyeData, i_BlinkDuringStim, i_NoFixation, criticalFixationData] = AnalyseEyeTrackerStuff2(eyeDataAllBlocks)

%Initialize output variables
nTrialsTotal = size(eyeDataAllBlocks.trial,1);
AllTrialNrs = 1:nTrialsTotal;

i_NoEyeData = false(nTrialsTotal,1);
i_BlinkDuringStim = false(nTrialsTotal,1);   
i_NoFixation = false(nTrialsTotal,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find trials in which data was missing or 'blinks/saccades' were detected during the critical period %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StimulusPeriods = [1 1.017; 1.517 1.534];   
FixationPeriod = [0.517 1.534];

for i=AllTrialNrs
    
    %Mark periods of no signal (zeros for longer than 1 second) as NaNs
    if ~isempty(eyeDataAllBlocks.events{i,1}.blinks.duration)
        idx_noData = find(eyeDataAllBlocks.events{i,1}.blinks.duration > 1000);
        for j=1:numel(idx_noData)
            timeStampsTemp = eyeDataAllBlocks.trial{i,1}(1,:);
            [~,idx_start] = min(abs(timeStampsTemp-eyeDataAllBlocks.events{i,1}.blinks.start_time(idx_noData(j))));
            [~,idx_end] = min(abs(timeStampsTemp-eyeDataAllBlocks.events{i,1}.blinks.end_time(idx_noData(j))));
            eyeDataAllBlocks.trial{i,1}(2:4,idx_start:idx_end) = NaN;
        end
    end
    
    %Mark the trials in which data is missing during the fixation period
    [~,idx_start] = min(abs(eyeDataAllBlocks.time{i,1}-FixationPeriod(1)));
    [~,idx_end] = min(abs(eyeDataAllBlocks.time{i,1}-FixationPeriod(2)));
    criticalSamples = isnan(eyeDataAllBlocks.trial{i,1}(2:4,idx_start:idx_end));
    if any(criticalSamples(:))
        i_NoEyeData(i) = true;
    end
    
    %Mark the trials in which there was a blink or saccade during stimulus presentation
    [~,idx_start1] = min(abs(eyeDataAllBlocks.time{i,1}-StimulusPeriods(1,1)));
    [~,idx_end1] = min(abs(eyeDataAllBlocks.time{i,1}-StimulusPeriods(1,2)));
    [~,idx_start2] = min(abs(eyeDataAllBlocks.time{i,1}-StimulusPeriods(2,1)));
    [~,idx_end2] = min(abs(eyeDataAllBlocks.time{i,1}-StimulusPeriods(2,2)));
    criticalTimeStamps = eyeDataAllBlocks.trial{i,1}(1,[idx_start1:idx_end1 idx_start2:idx_end2]);
    if ~isempty(eyeDataAllBlocks.events{i,1}.blinks.duration)
        criticalTimeStampsMatrix = repmat(criticalTimeStamps,[numel(eyeDataAllBlocks.events{i,1}.blinks.duration) 1]);
        blinkStartMatrix = repmat(eyeDataAllBlocks.events{i,1}.blinks.start_time,[1 numel(criticalTimeStamps)]);
        blinkEndMatrix = repmat(eyeDataAllBlocks.events{i,1}.blinks.end_time,[1 numel(criticalTimeStamps)]);
        i_blinkDuringStimTemp = any((blinkStartMatrix <= criticalTimeStampsMatrix)' & (criticalTimeStampsMatrix <= blinkEndMatrix)')';
        if sum(i_blinkDuringStimTemp) > 0
            i_BlinkDuringStim(i) = true;
        end
    end %The above was for blinks, now do the same for saccades (saccades often 'surround' the blinks)
    if ~isempty(eyeDataAllBlocks.events{i,1}.saccades.duration)
        criticalTimeStampsMatrix = repmat(criticalTimeStamps,[numel(eyeDataAllBlocks.events{i,1}.saccades.duration) 1]);
        blinkStartMatrix = repmat(eyeDataAllBlocks.events{i,1}.saccades.start_time,[1 numel(criticalTimeStamps)]);
        blinkEndMatrix = repmat(eyeDataAllBlocks.events{i,1}.saccades.end_time,[1 numel(criticalTimeStamps)]);
        i_saccadeDuringStim = any((blinkStartMatrix <= criticalTimeStampsMatrix)' & (criticalTimeStampsMatrix <= blinkEndMatrix)')';
        if sum(i_saccadeDuringStim) > 0
            i_BlinkDuringStim(i) = true;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find trials in which fixation was outside the bounds during the fixation period %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

width_of_1pix = eyeDataAllBlocks.widthOfScreen/eyeDataAllBlocks.screenResolution(1,1);                    %in meters 
height_of_1pix = eyeDataAllBlocks.heightOfScreen/eyeDataAllBlocks.screenResolution(1,2);                  %in meters 

criticalFixationData = cell(nTrialsTotal,4);
for i=AllTrialNrs
    
    if ~isempty(eyeDataAllBlocks.events{i,1}.fixation.duration)
    
        %Find median value of this trial
        nFixationPeriodsTemp = numel(eyeDataAllBlocks.events{i,1}.fixation.duration);
        x_avg_weight_temp = nan(nFixationPeriodsTemp,2);
        y_avg_weight_temp = nan(nFixationPeriodsTemp,2);
        for k=1:nFixationPeriodsTemp
            timeStampsTemp = eyeDataAllBlocks.trial{i,1}(1,:);
            [~,idx_start_temp] = min(abs(timeStampsTemp-eyeDataAllBlocks.events{i,1}.fixation.start_time(k)));
            [~,idx_end_temp] = min(abs(timeStampsTemp-eyeDataAllBlocks.events{i,1}.fixation.end_time(k)));
            x_relevant = eyeDataAllBlocks.trial{i,1}(2,idx_start_temp:idx_end_temp);
            y_relevant = eyeDataAllBlocks.trial{i,1}(3,idx_start_temp:idx_end_temp);
            x_relevant = x_relevant((~isnan(x_relevant)) & (x_relevant ~= 0));
            y_relevant = y_relevant((~isnan(y_relevant)) & (y_relevant ~= 0));
            if ~isempty(x_relevant)
                x_avg_weight_temp(k,:) = [median(x_relevant) numel(x_relevant)];
            end
            if ~isempty(y_relevant)
                y_avg_weight_temp(k,:) = [median(y_relevant) numel(y_relevant)];
            end
        end
        x_totalSamples = nansum(x_avg_weight_temp(:,2));
        x_avg_weight_temp(:,2) = x_avg_weight_temp(:,2) ./ x_totalSamples;
        y_totalSamples = nansum(y_avg_weight_temp(:,2));
        y_avg_weight_temp(:,2) = y_avg_weight_temp(:,2) ./ y_totalSamples;
        xRefTemp = nansum(x_avg_weight_temp(:,1).*x_avg_weight_temp(:,2));
        yRefTemp = nansum(y_avg_weight_temp(:,1).*y_avg_weight_temp(:,2));
        
        %Use the above median as a reference point for the samples in the critical fixation period of this trial    
        [~,idx_start_critical] = min(abs(eyeDataAllBlocks.time{i,1}-FixationPeriod(1)));
        [~,idx_end_critical] = min(abs(eyeDataAllBlocks.time{i,1}-FixationPeriod(2)));
        criticalSampleIdx = idx_start_critical:idx_end_critical;
        for j=1:nFixationPeriodsTemp   
            timeStampsTemp = eyeDataAllBlocks.trial{i,1}(1,:);
            [~,idx_start_temp] = min(abs(timeStampsTemp-eyeDataAllBlocks.events{i,1}.fixation.start_time(j)));
            [~,idx_end_temp] = min(abs(timeStampsTemp-eyeDataAllBlocks.events{i,1}.fixation.end_time(j)));
            idx_critical = criticalSampleIdx((idx_start_temp <= criticalSampleIdx) & (criticalSampleIdx <= idx_end_temp));
            criticalFixationData{i,1} = [criticalFixationData{i,1} idx_critical];
            criticalFixationData{i,2} = [criticalFixationData{i,2} eyeDataAllBlocks.time{i,1}(1,idx_critical)];
            
            x_relevant = eyeDataAllBlocks.trial{i,1}(2,idx_critical);
            y_relevant = eyeDataAllBlocks.trial{i,1}(3,idx_critical);
            
            x_relevant_pix = x_relevant-xRefTemp;
            y_relevant_pix = y_relevant-yRefTemp;
            
            %Convert the x and y gaze to degrees away from the fixation points
            x_relevant_angle = atand(x_relevant_pix .* width_of_1pix ./ eyeDataAllBlocks.distToScreen);
            y_relevant_angle = -1*atand(y_relevant_pix .* height_of_1pix ./ eyeDataAllBlocks.distToScreen);      %Note that we flip the axis here: positive angles are eye gazes towards the top
            criticalFixationData{i,3} = [criticalFixationData{i,3} x_relevant_angle];
            criticalFixationData{i,4} = [criticalFixationData{i,4} y_relevant_angle];
            
            %Mark the trials in which the participant did not fixate within 3 degrees of the reference point
            if any(abs(x_relevant_angle) > 3) || any(abs(y_relevant_angle) > 3)
                i_NoFixation(i) = true;
            end
        end
    end
end

return %[EOF]
