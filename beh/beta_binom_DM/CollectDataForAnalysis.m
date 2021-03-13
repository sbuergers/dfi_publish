function CollectedDataVariables = CollectDataForAnalysis(ResultsAllBlocks)

%Initialize output collection structure
CollectedDataVariables.dummy = [];

%Order the blocks based on the given blockNr in stimulus_info
nBlocks = size(ResultsAllBlocks,1);
Data = ResultsAllBlocks;
[~, indices_order] = sort(cell2mat(ResultsAllBlocks(:,6)));
for i=1:nBlocks
    for j=1:6
        Data{i,j} = ResultsAllBlocks{indices_order(i),j};
    end
end
clear ResultsAllBlocks indices_order i j

%Collect positional info on all trials
nTrialsPerBlock_planned = zeros(1,nBlocks);             %Collect the planned nr of trials per block
nTrialsPerBlock = zeros(1,nBlocks);                     %Collect the actually completed nr of trials per block
BlockNrs = [];                                          %Temporal position of the block in the whole experiment (in ResultsAllBlocks)
TrialNrs = [];                                          %Temporal position of the trial in the block (in behavioral_results)
for i=1:nBlocks
    stimulus_info = Data{i,1};
    nTrialsPerBlock_planned(i) = stimulus_info.nTrials;
    behavioral_results = Data{i,3};
    for j=1:length(behavioral_results)
        if ~isempty(behavioral_results{j,1})            %if the trial was completed
            nTrialsPerBlock(i) = nTrialsPerBlock(i)+1;
            BlockNrs = [BlockNrs; i];
            TrialNrs = [TrialNrs; j];
        end
    end   
end
nTrialsTotal_planned = sum(nTrialsPerBlock_planned);    %nTrials_planned of all blocks together
nTrialsTotal = sum(nTrialsPerBlock);                    %nTrials of all blocks together

%Collect some other info from the latest block's stimulus_info
experimentNr = stimulus_info.experimentNr;
experimentName = stimulus_info.experimentName;
Subj_nr = stimulus_info.Subj_nr;
taskName = stimulus_info.taskName;
eyeTrackerPresent = stimulus_info.eyeTrackerPresent;
distToScreen = stimulus_info.cpu.dist_to_screen;        %These sizes are used to compute visual angles of the eyeTracker results
widthOfScreen = stimulus_info.cpu.width_of_screen;
heightOfScreen = stimulus_info.cpu.height_of_screen;
clear stimulus_info behavioral_results i j

%Find conditions info on all trials
StimType               = nan(nTrialsTotal,1);           %Trial types (1-5)
ProbeLocation          = nan(nTrialsTotal,1);           %Locations in degrees visual angle
AuditoryLocation       = nan(nTrialsTotal,1);           %NaN = 'N/A'
VisualLocation         = nan(nTrialsTotal,1);           %NaN = 'N/A'
VisReliability         = nan(nTrialsTotal,1);           %The horizontal SD of the blob/cloud of dots
VisVariableRel         = nan(nTrialsTotal,1);           %Is the visual stimulus from fixed or variable reliability
LeadInTiming           = nan(nTrialsTotal,1);           %The lead-in time in milliseconds
for i=1:nTrialsTotal
    stim_spec = Data{BlockNrs(i),2};
    StimType(i)         = stim_spec{TrialNrs(i),1}.stim_type;
    ProbeLocation(i)    = stim_spec{TrialNrs(i),1}.probe_loc;
    AuditoryLocation(i) = stim_spec{TrialNrs(i),1}.aud_loc;
    VisualLocation(i)   = stim_spec{TrialNrs(i),1}.vis_loc;
    VisReliability(i)   = stim_spec{TrialNrs(i),1}.vis_horSD;
    VisVariableRel(i)   = stim_spec{TrialNrs(i),1}.variable_rel;
    LeadInTiming(i)     = stim_spec{TrialNrs(i),1}.timing_lead_in;
end
clear i

%Collect the responses
LocResponse     = nan(nTrialsTotal,1);            %-1 for Left, 1 for Right
LocResponseTime = nan(nTrialsTotal,1);            %in milliseconds
for i=1:nTrialsTotal
    behavioral_results = Data{BlockNrs(i),3}{TrialNrs(i),1};
    if isfield(behavioral_results,'LocResponse')
        if ~isempty(behavioral_results.LocResponse)             %if empty, then there was no answer, and LocResponse(i) remains NaN
            LocResponse(i) = behavioral_results.LocResponse;    
        end
    end
    if isfield(behavioral_results,'LocResponseTime')
        if ~isempty(behavioral_results.LocResponseTime)
            LocResponseTime(i) = behavioral_results.LocResponseTime;
        end
    end
end

%Collect some timing measurements
TrialDurations            = nan(nTrialsTotal,1);                            %in seconds
SaveDurations             = nan(nTrialsTotal,1);                            %in seconds
for i=1:nTrialsTotal
    behavioral_results = Data{BlockNrs(i),3}{TrialNrs(i),1};
    TrialDurations(i) = behavioral_results.totalTrialDuration;
    SaveDurations(i) = behavioral_results.timing_trialSaveTime;  
end
clear behavioral_results i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Collect some recurrent info and settings %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find the unanswered questions
nUnansweredTrials = sum(isnan(LocResponse));
    
%Find the unimodal trials
i_unimodal = isnan(VisualLocation) ~= isnan(AuditoryLocation);

%Find the AVdisparities that were used (note that we count congruent trials too: AVdelta = 0)     
AVdisparities = unique(VisualLocation(~i_unimodal) - AuditoryLocation(~i_unimodal));                      %visual is moved by +0.5 deltaAV, auditory is moved by -0.5 deltaAV

%Collect all locations
Locations = unique(ProbeLocation(~isnan(ProbeLocation)));
nLocations = numel(Locations);
for i=1:nLocations
    LocationLabels{i} = num2str(Locations(i));
end

%Collect the visual reliability levels
Reliabilities = unique(VisReliability(~isnan(VisReliability)));
nReliabilities = numel(Reliabilities);
for i=1:nReliabilities
    ReliabilityLabels{i} = num2str(Reliabilities(i));
end

%Collect the available stimulus types
StimTypes = unique(StimType);
nStimTypes = numel(StimTypes);
StimTypesNames = {'Aud only', 'Vis only', 'AV congruent', 'AV incongr VL AR', 'AV incongr VR AL'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Collect all variables into one large structure that serves as the output of this function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%Collect all variables in struct 'CollectedDataVariables'
clear CollectedDataVariables;                                           %clear the dummy!
varNames = who;
for i=1:length(varNames)
    eval(['CollectedDataVariables.' varNames{i} '=' varNames{i} ';']);
end

end %[EOF]

