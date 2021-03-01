function sys = setup_monitor(setup)
% setup = setup_monitor(setup)
%
% setup.name specifies which predefined screen setup to use ('office',
% 'lab', 'mac', 'macSteffen').
%
% The function returns setup with the additional field 'sys', containing 
% monitor specifications such as
%         sys.res = [1920 1080];
%         sys.refresh = 60;
%         sys.dist = 50; 
%         sys.width = 48; 
%         sys.fov
%         sys.ppd       
%
% ----------------------------------------
% adapted from Agoston Mikhalic, 
% last updated, April 2015
%

% Screen native resolutions and refresh rates, monitor width and distance params
switch setup.name
    case 'office'
        sys.res = [1920 1080];
        sys.refresh = 60;
        sys.dist = 50; % roughly
        sys.width = 48; % IIYAMA
    case 'lab'
        sys.res = [1920 1080];
        sys.refresh = 60;
        sys.dist = 50; % 50
        sys.width = 48; % IIYAMA
    case 'mac'
        sys.res = [1440 900];
        sys.refresh = [];
        sys.dist = 50; % roughly
        sys.width = 33; % roughly
    case 'macSteffen'
        sys.res = [1680 1050];
        sys.refresh = [];
        sys.dist = 50; % roughly
        sys.width = 33; % roughly
        sys.height = 21.5; % roughly
end

% Field of view (in degree) and number of pixels per cm
sys.fov  = atand(sys.width / 2 / sys.dist) * 2;
sys.ppd  = sys.res(1) / sys.fov; % inaccurate!

% eof
%           sys.width (in cm)
% ----S-----------------------
%      .  opp  |
%       .      |
%        .     |
%         .    |
%          .   |sys.dist (in cm)
%           .  |adj
%   vis angle. |
%     (degree).|
%              |
%          headrest



