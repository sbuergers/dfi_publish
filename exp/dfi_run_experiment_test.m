function dfi_run_experiment
% dfi_run_experiment
%
% GENERAL:
%   Command line short-cut for an audiovisual presentation paradigm that
%   presents flashes of light and auditory clicks, which might induce a 
%   double flash illusion or fusion effects.
%
%
% PARADIGM (yesno - EEG):                                           
%  |----------------------------|...............|----------X---------------|
% strl    ITI(alpha)           AV      <fe>     V          R      slack   etrl
%                              AV      <dfi>   A(V)
%
%         ITI     = unif(0.6, 1.6)
%         StimWin = {0.0167:0.192} = 0.0771 
%         el      < 0.02
%         Resp    = 1.5 (constant)
%         Slack   = 1.5         
%         AVERAGE TRIAL DURATION <= 4.2 s
%         CONDITIONS: F   (2)
%                     S   (2)
%                     FS  (4)
%         ncond * ntrls * nsoa * trldur = 8*100*8*3.5 = 373.4 minutes 
%         which gives 75 minutes of pure testing for 5 sessions total,
%         adding 20 minutes for breaks etc. gives 95 minutes + EEG
%         preparation, which is roughly 75 minutes and conservatively 85
%         minutes, giving 180 minutes or 3h / session
%         Though I'm not sure this can be achieved throughout, so to be
%         safe I should start out with 3:30h or even 4h.
%
%
% PARADIGM (2IFC - behavioral only):                                           
%  |---------------|......|---------------|----------------X---------------|
% strl    ITI     Va1    Va2    2IFCI    Vb1               R      slack   etrl
%                        
%
%         ITI     = unif(0.6, 1.6)
%         2IFCI   = 0.8
%         StimWin = 0.0771 + 0.0167 = 0.0938
%         Resp    = 1.5 (exactly)
%         Slack   = 1.5         
%         AVERAGE TRIAL DURATION <= 5.2 s
%         DESIRED NUMBER OF TRIALS PER CONDITION: min. 40, so try 45
%         CONDITIONS: 2F   - 1F   and vice versa
%                     2F1S - 1F1S and vice versa 
%                     2F2S - 1F2S and vice versa 
%                       2S -   1S and vice versa
%         ncond * ntrls * nsoa * trldur = 4*50*8*5.2 < 139 minutes 
%         Plus some spare time for setting up, breaks etc, say 10 minutes
%         per run, so 160 minutes over two sessions (70 min each, +- 5 min)
%
%
%                               
% DOUBLE FLASH ILLUSION:
%   Physical stimuli often emmit both light and sound, and as human
%   perception is inherently a noisy representation of the environment it
%   can be highly advantageous to use information from multiple senses to
%   determine how to act. In many cases this results in perceptual
%   integration of sound and vision into a more reliable internal
%   repjresentation of a physical stimulus. However, this can also lead to 
%   perceptual distortions and illusions. In the double flash illusion (DFI) 
%   a visual stimulus is paired with an auditory stimulus after which a second 
%   auditory stimulus is presented in close succession. In some instances,
%   this second tone is able to elicit the percept of a concurrent visual 
%   stimulus. This is sometimes called fission. Cecere et al. (2014) have
%   shown that the DFI window of integration is directly correlated with
%   individual alpha peak frequency, with longer wavelengths leading to
%   larger windows of integration. 
%
% FUSION EFFECT:
%   In contrast to the DFI the fusion effect (FE) emerges from unisensory
%   failure to separate two distinct stimuli, like two successive pulses of
%   light. At a certain inter stimulus interval, people will tend to
%   sometimes integrate two flashes into the percept of a single flash
%   whereas on other trials two stimuli will be perceived. Interestingly,
%   the window for the FE appears to be similar to the window of the DFI
%   and has been associated with EEG/MEG recordings of occipital alpha
%   frequencies (Lange et al., 2013). Also, the original FE was described
%   by Andersen 2004, and several times replicated (and not replicated)
%   since. In this effect a single auditory stimulus can increase the
%   likelihood of perceiving a single flash as opposed to two flashes. 
%
% RESEARCH QUESTION:
%   It has been found that alpha-oscillations influence perception in
%   various ways (Hanslmayr et al., 2011). Alpha power is associated with 
%   increased internal processing and reduced external processing, leading 
%   to reduced perceptual sensitivity. Similarly, alpha up-states are 
%   associated with external, whereas down-states are associated with 
%   internal processing. Finally, alpha phase-coherence between cortical 
%   sites is associated with augmented internal processing as well. 

%
%   Cecere et al. (2014) showed that in a double flash illusion paradigm,
%   individual alpha frequency was associated with the window of
%   integration and the perception an illusory flash, with broader wavelengths
%   being associated with larger windows of integration. They furthermore
%   found that individual alpha power was inversely related to perceiving
%   the illusion (high power ~ low illusion probability). 
%   This is in line with the results of Lange et al. (2013), where the 
%   DFI is associated with low pre-stim alpha power and the FE with high 
%   alpha power. The authors interpreted this as meaning that low alpha
%   power in the visual cortex signifies enhanced excitability and
%   therefore a higher likelihood of perceiving two flashes instead of one
%   in the case of the FE paradigm, but also to perceive two flashes 
%   instead of one in the DFI, as internal or cross-modal input, similar to
%   external input is more likely to excite visual neural assemblies. 
%
%   We want to know whether the effect found by Cecere et al. (2014) holds
%   within subjects and whether this is specific to the DFI or if we can 
%   also find this for the fusion effect. Finally, it might be of interest 
%   to use TMS to look for entrainment effects of alpha frequencies and a 
%   direct causal relationship. 
%   
% USAGE:
% Parameters to be specified within code. 
%
% NOTES:
% You might want to try to run this from the command line and use the 
% matlab -nojvm flag to disable unnecessary graphics operations for more
% precise stimulus timing and system respondence. This is untested so far
% though, and possibly unnecessary.
%
%
% ========================================================================
% EXTRA CONDITIONS
%
% MULTIMODAL CONDITION ATTEND AUDITORY:
% I can run the multimodal attend auditory condition simply by specifying
% multimodal as condition and decreasing the SOA to the lowest SOA
% possible (plus changing the number of trials and blocks).
%
% PERFECT PSE RUNS:
% Select only the SOA giving the best approximation to the PSE (roughly 50/50
% responses). Select trial number (e.g. 300 for multimodal and 
%
% Calls the functions:
%   dfi_setup_dir.m
%   dfi_setup_monitor.m
%   dfi_getaudio.m
%   dfi_generate_data.m
%   dfi_present_stimuli.m
%
% ----------------------------------------
% adapted from Agoston Mihalik, 
% last updated, September 2015
%

% Initial clean up 
if feature('IsDebugMode'), dbquit; end;
clc
close all
clear functions
clearvars
commandwindow;

% add java garbage collector and check working directory
if exist('dfi_run_experiment.m', 'file')
    addpath(genpath(pwd))
    javaaddpath(fullfile(pwd, 'jheapcl', 'MatlabGarbageCollector.jar'));
else
    error('ERROR: Please go to the experiment folder "dfi_exp_Steffen".')
end


%% -------------------------- start settings -------------------------- %%
s.practicefeedback = 0; % leave this off, otherwise feedback is given after each trial

% Debug mode 
dbstop if error  
debug = 0;    
s.stim.plot = debug; 

% Random number generator seed (shuffle = random seed)
% if ~verLessThan('matlab', '8.3'), rng('shuffle');  end 



%% 0.)  Set Feedback 
% Note that you always have to also check physical timing, i.e. with
% photodiode and microphone, and not rely solely on PTB flip times.
s.fb.flipdur       = 0;     % fb only works in yesno paradigm.
s.fb.flipdiagbig   = 0;
s.fb.flipdiagsmall = 1;
s.fb.plotres       = 1;



%% 1.)  Input prompt and path 
prompt  = {'Setup name', 'Subject id', 'Condition', 'Pilot', 'Session id', 'Run id', 'Include Practice', 'Paradigm', 'Response key order', 'Staircase'};
title   = 'Input parameters';
lines   = [1 35]; % allow 35 chars in each input box
default = {'PiaCRT', '999', 'all', '1', '1', '1', '1', 'yesno', 'AA', 'no'}; % (A = left 1, right 2, B = left 2, right 1) 
params  = inputdlg(prompt, title, lines, default);
if ~isempty(params) % ok 
    [s.name, s.subj.id, s.cond, s.paradigm, s.keyswitch, s.stair.incl] = params{[1 2 3 8 9 10]};
    param_strings = cellfun(@str2num, params([4 5 6 7]), 'UniformOutput', false);
    [s.subj.pilot, s.subj.session, s.runid, s.practice] = param_strings{:};
    if strcmp(s.keyswitch, 'AA') && ~strcmp(s.cond, 'multiAttAud') && ~strcmp(s.paradigm, '2AFC')
        error('All conditions except multiAttAud require counterbalancing of response keys, specify AB or BA (A = left 1, right 2, B = left 2, right 1) depending on previous run')
    end
else % cancel 
    sprintf('Run has been aborted...\n')
    return
end
% add psychtoolbox path for pia
if strcmp(s.name, 'PiaCRT') || strcmp(s.name, 'SteveCRT')
    addpath /Applications/Psychtoolbox;
end
s.mypath = dfi_setup_dir(s); 
if strcmp(s.paradigm, '2IFC'), s.paradigm = '2AFC'; end;

% keep diary of command line history to see warnings etc.
time = datestr(now, 'HH:MM');
diary(sprintf('cmdline_history_subj%s_sess%g_run%g_%s.mat', ...
               s.subj.id, s.subj.session, s.runid, [time(1:2),'h',time(4:5),'m']))
diary on



%% 2.)  Display
s.disp = dfi_setup_monitor(s);
s.disp.incolor = 0; % gray-scale or red over green?
% Text font and size
s.text.font = 'Arial';
s.text.size = 16;



%% 3.)  Response devices (Keyboard and mouse) 
s.respdev = 'keyboard'; % keyboard, mouse

% Keyboard 
[kbids, kbnames] = GetKeyboardIndices;
if numel(kbids) < 2, warning('There is no external keyboard connected!'); end;
s.kb.subj = kbids(strcmp(kbnames,'Wired Keyboard 400'));                 % participant
s.kb.op   = kbids(strcmp(kbnames,'Apple Internal Keyboard / Trackpad')); % operator
% resp and quit keys
s.key.resp  = {'F' 'J'}; % 1 and 2 is for numpad, the other 1 and 2 are '1!' and '2@' respectively
if strcmp(s.paradigm, 'YN_threshold')
    s.key.resp  = {'A' 'S' 'D' 'F' 'J' 'K' 'L' ';:'};
end
s.key.instruction = {'LEFT' 'RIGHT'}; 
s.key.quit        = {'q'};

% Mouse (always returns the unified state of all mice, ids don't exist)



%% 4.)  Eyetracker
s.eyelink.enable    = 0;
s.eyelink.dummy     = 0;
s.edfFlag           = '-e';    % -e only convert event samples to asc format (edf still contains everything)
s.eyelink.fixdeg    = 5;       % move up until fixdeg degrees away from absoluate central position (512, 384)
s.eyelink.driftdeg  = 1;       % allowed discrepancy from center during drift check
s.eyelink.fb        = 'block'; % 'trial' or 'block' (or any other string for no fb)



%% 5.)  Stimuli
if strcmp(s.name, 'PiaCRT')
    s.stim.vposdelay = -0.0025; %-0.0029327;  
elseif strcmp(s.name, 'SteveCRT')
    s.stim.vposdelay = -0.001;  
end
s.stim.vdur  = 0.008;   %   12 ms (Cecere), 17 ms (Shams), 1.3 ms (Innes-Brown 2009) CRT:~2ms
s.stim.adur  = 0.002;   %    7 ms (Cecere)                       
s.stim.afreq = 3500;    % 3500 Hz (Cecere)                      
refr         = 1/s.disp.refr;   
if strcmp(s.paradigm, 'yesno') 
    % s.stim.SOA     = [3, 5, 6, 7, 9, 13, 25]*refr; % pilot 05 and 06 of May
    s.stim.SOA     = [3, 5, 6, 7, 9, 13, 19, 27]*refr; % yesno
elseif strcmp(s.paradigm, '2AFC') 
    %s.stim.SOA     = [3, 4, 6, 7, 9, 13, 30]*refr; % 2IFC
    %s.stim.SOA     = [3, 4, 6, 7, 9, 13, 30]*refr; % 2IFC
    %s.stim.SOA     = [3, 5, 6, 7, 9, 13, 27]*refr; % 2IFC // participant 701, ie May
    s.stim.SOA     = [3, 5, 6, 7, 9, 13, 19, 27]*refr; % yesno
elseif strcmp(s.paradigm, 'YN_threshold')
    % assign custom SOAs for illusion windows
    prompt  = {'FF (in flips!)', 'Fusion illusion (in flips!)', 'Fission illusion (in flips!)', 'FSFS (in flips!)'};
    title   = 'SOAs (subject specific)';
    lines   = [1 35]; % allow 35 chars in each input box
    default = {'', '', '', ''};
    params  = inputdlg(prompt, title, lines, default);
    if ~isempty(params) % ok
        soas_YN_threshold = cellfun(@str2num, params([1 2 3 4]), 'UniformOutput', false);
        s.stim.SOA_FF  = soas_YN_threshold{1}*refr;
        s.stim.SOA_Fus = soas_YN_threshold{2}*refr;
        s.stim.SOA_Fis = soas_YN_threshold{3}*refr;
        s.stim.SOA_2av = soas_YN_threshold{4}*refr;
    end
    s.stim.SOA = 1; % for calculation of number of trials later on
end
s.stim.SOA2AFC = 0.8; % 0.9;
s.stim.vloc    = [15 0];  % [x y]                                    % \yesno\
                                                                     % id  VA
% nopbook2:0.005878 + 0.002358; % for vloc = [15 0]                  % 1 - 00
if strcmp(s.paradigm, 'yesno')                                       % 2 - 10  = V1
    if strcmp(s.cond, 'multimodal')                                  % 3 - 20  = V1V2
        s.stim.trlid = [5 6 8 9];                                    % 4 - 01  = A1
    elseif strcmp(s.cond, 'visual')                                  % 5 - 11  = A1V1
        s.stim.trlid = [2 3 2 3];                                    % 6 - 21  = V2A1
    elseif strcmp(s.cond, 'audio')                                   % 7 - 02  = A1A2
        s.stim.trlid = [4 7 4 7];                                    % 8 - 12  = V1A2
    elseif strcmp(s.cond, 'all')                                     % 9 - 22  = V2A2 
        s.stim.trlid = [2 3 5 6 8 9 2 3 5 6 8 9]; 
        % nil trials are added in dfi_generate_data
    elseif strcmp(s.cond, 'multiAttAud')
        s.stim.trlid = [5 6 8 9];
    end                                                         
elseif strcmp(s.paradigm, '2AFC') || strcmp(s.paradigm, '2AFCsoundTest') % \2AFC\
    if strcmp(s.cond, 'multimodal')                             % 2 -   V1-V1V2  =  V
        s.stim.trlid = [5 6 8 9 ];                              % 3 - V1V2-V1    =  V
    elseif strcmp(s.cond, 'visual')                             % 4 -   A1-A1A2  =  A
        s.stim.trlid = [2 3 2 3];                               % 5 - V1A1-V2A1  =  Fus
    elseif strcmp(s.cond, 'audio')                              % 6 - V2A1-V1A1  =  Fus
        s.stim.trlid = [4 7 4 7];                               % 7 - A1A2-A1    =  A
    elseif strcmp(s.cond, 'all')                                % 8 - V1A2-V2A2  =  Fis
        s.stim.trlid = [2 3 5 6 8 9 2 3 5 6 8 9];               % 9 - V2A2-V1A2  =  Fis
    elseif strcmp(s.cond, 'multiAttAud')
        s.stim.trlid = [5 6 8 9];
    end                                                         
elseif strcmp(s.paradigm, 'YN_threshold')      % Only present illusion trials 
    s.stim.trlid = [2 3 5 6 8 9];
end

                                                                 

%% 6.)  Trial / Design                
% PARADIGM:            ready           soa                             
%  |--------------------|.......|...............|----------X---------------|
% strl    ITI(alpha)            V      <fe>     V          R      resw    etrl
%                              AV      <dfi>   A(V)  
% experimental parameters
% Note that atm you have to choose very specific block and trial settings
% in order for everything to be balanced. So Nav = 2Nv, Na = 0. 
s.ITI        = [0.25 0.25];   % unif(m, sd) measure alpha here
s.reswin     = 1.5 + 0.6;           % max response window
if strcmp(s.paradigm, '2AFC')
    s.reswin = 1.5; 
    if strcmp(s.respdev, 'mouse') % give extra second for confidence judgment
        s.reswin = 2.5;
    end
end
if strcmp(s.paradigm, '2AFC')
    s.slack = 1;
else
    s.slack = 1.2; %1.5;      % let response related activity taper out
end
s.pseudoran  = 1;             % if 1, prohibit randomly occuring streaks
s.nallowrep  = 2;             % of one trial type with soa x
s.breakint   = 3;             % forced break after every n blocks
s.breakdur   = 20;            % forced break duration in s
s.randwblocks= 1;             % randomize blocks, or all trl types within blocks (=1)?
s.ntrls      = numel(s.stim.SOA)*numel(s.stim.trlid) + 1; % per block (plus one nil trial)
    s.blocks.N   =  3; % nblocks                % ntrls    soa    conds
    s.blocks.Nv  =  4;                          % 100   *   7   *   7   =   4900
    s.blocks.Nav =  6;                          % Session 1: N = 15
s.blocks.Na  = 0;
if sum([s.blocks.Nv, s.blocks.Na, s.blocks.Nav]) ~= s.blocks.N && ~s.randwblocks==1
    error('Block number of individual modalities does not add up to total block number.')
end
% practice parameters
s.pract.ntrls      = 24;
s.pract.blocks.N   = 2; % number of visual or multimodal practice blocks
s.pract.blocks.Nv  = 1;
s.pract.blocks.Nav = 1;
s.pract.blocks.Na  = 0;

% make sure you don't have to change everything by hand all the time
if strcmp(s.paradigm, '2AFC')
    s.ntrls      = numel(s.stim.SOA)*numel(s.stim.trlid); % per block
    s.blocks.N   =  12;  
    if strcmp(s.cond, 'visual')
        s.blocks.N = 5;
    end % Session 4: N = 6
%     s.blocks.Nav =  6;                       
end

% auditory only condition with shortest soa
if strcmp(s.cond, 'multiAttAud')
    s.stim.SOA = [3]*refr; % yesno
    s.ntrls = numel(s.stim.SOA)*numel(s.stim.trlid)*12;
%     s.ntrls = numel(s.stim.SOA)*numel(s.stim.trlid)*10;
    s.blocks.N = 4;
    s.nallowrep = 4;
end

% YN_threshold (optimal illusion window SOAs)
if strcmp(s.paradigm, 'YN_threshold')
    s.reswin     = 2.5 + 0.6; % 2.1 ...
    s.nallowrep  = 4;
    s.blocks.N   = 2;
    s.ntrls      = 12 * 6 + 1; % 2 nil trials per block
    s.stim.trlid = [2 3 5 6 8 9]; 
end



        
%% 7.)  Audio settings
% Sound sampling frequency (needs to be high to have accurate ITDs)
s.stim.aud.freq = 192000; % 44100
s.stim.aud.freq = 96000;  % 44100
% s.stim.aud.freq = 48000;  % 44100
s.stim.aud.vol  = 1;   %  // max speaker volume and 4 bar earphone volume
s.stim.aud.ramp = 1;     % Cecere/Shams didn't use ramp
s.stim.aud.rdur = 0.0005;
% ITD (interaural time difference)
s.stim.ITD.soundspeed = 343.2;   % speed of sound in dry air at 20 ?C in m/s.
s.stim.ITD.eardist    = 0.2;     % in m (approx.)
s.stim.ITD.timelag    = s.stim.ITD.eardist * (sind(s.stim.vloc(1)) + (s.stim.vloc(1) * pi / 180)) / ...
                        s.stim.ITD.soundspeed; % based on "psychology of hearing", untested
                 


%% 8.)  Video settings
% Visual stimulus parameters
s.stim.vis.type = 'blob'; % 'blob', 'circle'
s.stim.vis.jittercntr = 0;
s.stim.vis.contrbnds  = [0.15 0.3];
s.stim.vis.jitterpos  = 0; 
s.stim.vis.jitperc    = 0.04; % 1 = 100% of diameter in any direction
s.stim.vdiam     = [2 2];
if strcmp(s.stim.vis.type, 'blob') % gaussian blob
    s.stim.vdiam = [4 4];
end
if strcmp(s.paradigm, '2AFCsoundTest')
    s.ITI      = [0.5 0];  % (m, sd) measure alpha here
    s.stim.SOA2AFC = 0.20;
    s.reswin   = 0;        % max response window
    s.slack    = 0;        % let response related activity taper out
    s.blocks.N = 1;        % nblocks
    s.ntrls    = numel(s.stim.SOA)*numel(s.stim.trlid);    
    s.stim.aud.vol  = 1;
    s.stim.adur     = 0.002;
    s.stim.aud.ramp = 0;
    s.stim.vis.contr   = 0.5; % perc. lighter than bgcontr
    s.stim.vis.bgcontr = 0.5; % perc. of maximum white
else
    s.stim.vis.contr   = 0.5; % perc. lighter than bgcontr
    s.stim.vis.bgcontr = 0.5; % perc. of maximum white
%     s.stim.vis.visib = 0.8;
end



%% 9.)  LABJACK settings
if ~exist('s.lj','class')
    %open a U3 with verbose logging on 
    s.lj = labJack('deviceID',3,'verbose',true);
end
if strfind(s.lj.version,'FAILED')
    warning('Opening labjack failed! Unable to send triggers!');
    s.ljPresent = false;
else
    s.ljPresent = true;
end



%% 10.) Staircase settings
s.stair.nChains = 2;
s.stair.nTrials = 4;
% initialize staircase chains
if strcmp(s.stair.incl, 'yes')
    ud_general = PAL_AMUD_setupUD(...
        'up', 1, ...
        'down', 1, ...
        'stepSizeUp', 1/120, ...
        'stepSizeDown', 1/120, ...
        'stopRule', s.stair.nTrials);
    % starting values are estimated from previous YN data
    for ichains = 1:s.stair.nChains
        s.ud{ichains}.FF  = PAL_AMUD_setupUD(ud_general, 'startValue', s.stim.SOA_FF);
        s.ud{ichains}.Fus = PAL_AMUD_setupUD(ud_general, 'startValue', s.stim.SOA_Fus);
        s.ud{ichains}.Fis = PAL_AMUD_setupUD(ud_general, 'startValue', s.stim.SOA_Fis);
    end
    % track if all chains have converged
    s.stair.fullstop = 0;
end



%% 11.) EEG markers 
% labjack has 8 bits, i.e. 256 is max
if strcmp(s.paradigm, 'yesno') || strcmp(s.paradigm, 'YN_threshold')
    % stimuli A
    %        X1 = one flash
    %        X2 = two flashes , X denotes the type of stimulus
    s.trig.siti =  10;             % start of ITI (after response slack)        
    s.trig.sV   =  20;
    s.trig.sA   =  30;
    s.trig.sAV  =  40;             
    s.trig.sV2  =  50;
    s.trig.sA2  =  60;
    s.trig.sAV2 =  70;
    s.trig.nil  =  00;
elseif strcmp(s.paradigm, '2AFC') || strcmp(s.paradigm, '2AFCsoundTest')
    % stimuli A
    s.trig.siti   = 10;        
    s.trig.a.sV   = 20;   
    s.trig.a.sA   = 30;
    s.trig.a.sAV  = 40;             
    s.trig.a.sV2  = 50;
    s.trig.a.sA2  = 60;
    s.trig.a.sAV2 = 70;
    % Stimuli B
    s.trig.b.sV   = 20;
    s.trig.b.sA   = 30;
    s.trig.b.sAV  = 40;             
    s.trig.b.sV2  = 50;
    s.trig.b.sA2  = 60;
    s.trig.b.sAV2 = 70;
end
s.trig.rOne = 100;     % respose 'I saw one', OR 'Stimulus A had 2'
s.trig.rTwo = 200;     % respose 'I saw two',    'Stimulus B had 2'

s.trig.sBlock = 123;   % send at start of a block
s.trig.eBlock = 132;   % send at end of a block



%% 12.) Debugging specs, also syn test (uncomment for quick run)
%         % for testing a whole run very quickly - e.g. test instructions
%         s.ITI       = [0.1 0];%[0.6 0];    
%         s.stim.SOA2AFC = 0.20;
%         s.stim.SOA     = [2, 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]*refr;
%         s.stim.trlid   = [9, 9, 9, 9, 9];
%         s.reswin    = 0.1;           
%         s.slack     = 0.1;          
%         s.pseudoran = 0;
%         s.nallowrep = 3;
%         s.breakint  = 20;             
%         s.breakdur  = 0;            
%         s.stim.vis.contrbnds  = [0.5 1];
%         s.stim.vis.jittercntr = 0;
%         s.stim.vis.jitterpos  = 0; 
%         s.stim.vis.contr   = 0.5; % perc. lighter than bgcontr
%         s.stim.vis.bgcontr = 0.5; % perc. of maximum white
%         s.stim.aud.vol  = 1;   %  // max speaker volume and 4 bar earphone volume


% ---------------------------- end settings ---------------------------- %


%% 13.) Run paradigm

% setup
s = dfi_setupPTB(s, debug);        
s = dfi_setupKeyboard(s);  
s = dfi_setup_el(s); % eyetracker    

dfi_instructions('expStart', s) 
if ~strcmp(s.paradigm, '2AFC')
    dfi_instructions('expStart2', s) 
end
dfi_instructions('response_finger_position', s)
  
% practice run
if s.practice
    [s, dpractice] = dfi_present_stimuli(s, debug); 
end

% experiment
[s, d] = dfi_present_stimuli(s, debug);

dfi_instructions('endLastBlock', s)

% plot outcomes
if s.fb.plotres, dfi_plotResults(s, d); end; 

% clean up
cleanWorkspace = 0;
if exist('dpractice', 'var')
    dfi_cleanUp(s, d, cleanWorkspace, debug, dpractice);
else
    dfi_cleanUp(s, d, cleanWorkspace, debug);
end

% Clear matlab debug mode
dbclear if error

% eof










% Potentially useful code:
% 
% % Set the blend funciton for the screen - good for anti aliasing
% Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');


% GetSecsTest
 
% This test is meant for Microsoft Windows only!
%  
% Performs a reliability test of your systems timing hardware. This script
%  tries to find out if your systems clock works correctly, ie., if
%  GetSecs(), WaitSecs(), Screen('Flip') and the PsychPortAudio functions
%  for timed stimulus onset and clock queries will work correctly.









