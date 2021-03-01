function s = dfi_setup_el(s)
%s = dfi_setup_el(s)
% tries to connect to an EyeLink eyetracker. Performs the initialization of 
% eyelink, connection to the eyelink data file, as well as calibration and 
% drift correction. If this works out the eye tracker can be used during 
% stimulus presentation.
%
% INPUT:  s.eyelink.dummy    -    contains a logical, whether or not to use
%                                 dummy mode for initialization
%
% OUTPUT: s.el   -   structure containing specifications, and used by eyelink
%                    functions

% Sample data:
% Keyword     Data Type
% LEFT,       Sets the intended tracking eye (usually include both LEFT and...
% RIGHT       RIGHT)
% GAZE        includes screen gaze position data
% GAZERES     includes units-per-degree screen resolution at point of gaze
% HREF        head-referenced eye position data 
% HTARGET     target distance and X/Y position (EyeLink Remote only)
% PUPIL       raw pupil coordinates 
% AREA        pupil size data (diameter or area)
% BUTTON      buttons 1-8
% STATUS      state and change flags warning and error flags (not yet supported)
% INPUT       input port data lines (not yet supported)
%
% Event data:
% Keyword     Effect
% GAZE        includes display (gaze) position data.
% GAZERES     includes units-per-degree screen resolution (for start, end of event)
% HREF        includes head-referenced eye position
% AREA        includes pupil area or diameter
% VELOCITY    includes velocity of parsed position-type (average, peak, start and end)
% STATUS      includes warning and error flags, aggregated across event (not yet supported)
% FIXAVG      include ONLY averages in fixation end events, to reduce file size
% NOSTART     start events have no data other than timestamp
%
% GAZE,GAZERES,AREA,HREF,VELOCITY  - default: all useful data
% GAZE,GAZERES,AREA,FIXAVG,NOSTART - reduced data for fixations
% GAZE,AREA,FIXAVG,NOSTART         - minimal data
%   


% Do we even want to record eyetracking data?
if s.eyelink.enable == 0
    fprintf('\nNot recording eye-tracking data on this run.\n');
    s.el.online = 0;
    return;
end


% Settings
s.el = EyelinkInitDefaults(s.disp.win);

s.el.backgroundcolour        = s.disp.grey;
s.el.msgfontcolour           = s.disp.white;
s.el.imgtitlecolour          = s.disp.white;
s.el.calibrationtargetcolour = s.disp.white;
s.el.calibrationtargetsize   = 1;
s.el.calibrationtargetwidth  = 0.5;
s.el.targetbeep              = 0;
s.el.feedbackbeep            = 0;
s.el.displayCalResults       = 1;
s.el.eyeimagesize            = 50;  % percentage of screen
s.el.eye_used                = s.el.LEFT_EYE;
% s.el.allowlocaltrigger       = 0; % allow user to trigger him or herself
% s.el.allowlocalcontrol       = 0; % allow control from subject-computer

EyelinkUpdateDefaults(s.el);


% Initialization of the connection with the eyetracker.
s.el.online = EyelinkInit(s.eyelink.dummy, 1);
if ~s.el.online
    fprintf('Eyelink Init aborted.\n');
    cleanup;  % cleanup function
    return;
end


% edf link
if strcmp(s.paradigm, 'yesno')
    condition_initial = 'y';
elseif strcmp(s.paradigm, '2AFC')
    condition_initial = 'i';
else
    condition_initial = 't';
end
s.edfFile = sprintf('%s%s%d%d.edf',s.subj.id, condition_initial, s.subj.session, s.runid);
res = Eyelink('Openfile', s.edfFile);
Eyelink('Command', 'add_file_preamble_text = "DFI experiment recording of participant %s', s.subj.id);
if res~=0
    fprintf('Cannot create EDF file ''%s'' ', s.edfFile);
    cleanup;
    return;
end


% make sure we're still connected.
if Eyelink('IsConnected')~=1 && ~s.eyelink.dummy
    cleanup;
    return;
end


% This Command is crucial to map the gaze positions from the tracker to
% screen pixel positions to determine fixation
Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, s.disp.res(1)-1, s.disp.res(2)-1);
Eyelink('message','DISPLAY_COORDS %ld %ld %ld %ld',        0, 0, s.disp.res(1)-1, s.disp.res(2)-1);

% Use conservative online saccade detection (cognitive setting)
Eyelink('Command', 'recording_parse_type = GAZE');
Eyelink('Command', 'saccade_velocity_threshold = 30');
Eyelink('Command', 'saccade_acceleration_threshold = 9500');
Eyelink('Command', 'saccade_motion_threshold = 0.15');
Eyelink('Command', 'saccade_pursuit_fixup = 60');
Eyelink('Command', 'fixation_update_interval = 0');

% Other tracker configurations
Eyelink('Command', 'calibration_type = HV9');
Eyelink('Command', 'generate_default_targets = YES');
Eyelink('Command', 'enable_automatic_calibration = YES');
Eyelink('Command', 'automatic_calibration_pacing = 1000');
Eyelink('Command', 'binocular_enabled = NO');
Eyelink('Command', 'use_ellipse_fitter = NO');
Eyelink('Command', 'sample_rate = 2000');
Eyelink('Command', 'elcl_tt_power = %d', 3); % illumination, 1 = 100%, 2 = 75%, 3 = 50%


% set edf data
Eyelink('Command', 'file_event_filter = LEFT,FIXATION,SACCADE,BLINK,MESSAGE,INPUT');
Eyelink('Command', 'file_sample_data  = LEFT,GAZE,GAZERES,AREA,HREF,STATUS,INPUT');

% set link data (can be used to react to events online)
Eyelink('Command', 'link_event_filter = LEFT,FIXATION,SACCADE,BLINK,MESSAGE,FIXUPDATE,INPUT');
Eyelink('Command', 'link_sample_data  = LEFT,GAZE,GAZERES,AREA,STATUS,INPUT');


% Calibrate the eye tracker
EyelinkDoTrackerSetup(s.el);

% do a final check of calibration using driftcorrection
success=EyelinkDoDriftCorrection(s.el);
if success~=1
    cleanup;
    return;
end    
    
% Cleanup routine:
function cleanup
% Shutdown Eyelink:
    Eyelink('Shutdown');
    s.el.online = 0;
end


end % end fun
% eof

 
 
















