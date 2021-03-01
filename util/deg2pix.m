function pixel = deg2pix(degree, sys)
%deg2pix(degree, sys)
% Calculates pixel position from visual angle (degree) in a hemi-field (see
% picture). You can give it the screen 'width' or 'height' in a field of 
% sys or specify a window pointer field 'win' to let PTB estimate it.
% If 'win' is provided it is always used to estimate width and height, 
% irrespective of whether they are provided as well. Another field 'res'
% needs to include the screen resolution in pixels (e.g. [1680 1050]).
%
% Returns pixel(1)=horizontal and 
%         pixel(2)=vertical number of pixels
%
% Note that you have to call the function separately for different degrees
% for height and width.
%
%      sys.width (in cm)
% ----S<--pix-->----------------
%      .  opp  |
%       .      |                        FORMULA: opp = adj * tan(a)
%        .     |
%         .    |
%          .   |sys.dist (in cm)
%           .  |adj
%   vis angle. |
%     (a)     .|
%              |
%          headrest
%
% EXAMPLE:
%   sys.res      = [1680, 1050]; % screen resolution
%   sys.width    = 30; % cm
%   sys.height   = 20; % cm
%   sys.dist     = 80; % cm
%   eccentricity = [10 0]; % degrees
%   pix = deg2pix(eccentricity, sys)
%
% ----------------------------------------
% adapted from Agoston Mihalic, 
% last updated, April 2015
%

    % enough information?
    if ~isfield(sys, 'win')
        if ~isfield(sys, 'width') && ~isfield(sys, 'height')
            error('I do not know the size of your monitor. Please provide it or a window pointer for estimation.')
        else
            if isfield(sys, 'width'),  width  = sys.width;  end;
            if isfield(sys, 'height'), height = sys.height; end;
        end
    else
        [width, height]=Screen('DisplaySize', sys.win);
        width = width/10; % convert to cm 
        height= height/10;
    end
    % width
    if exist('width', 'var')
        ppcmw    = sys.res(1)/width;
        pixel(1) = tand(degree(1)) * sys.dist * ppcmw;
    end
    % height
    if exist('height', 'var')
        ppcmh    = sys.res(2)/height;  
        pixel(2) = tand(degree(2)) * sys.dist * ppcmh;
    end
end % end deg2pix
% eof
