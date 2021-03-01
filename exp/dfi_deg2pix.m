function pixel = dfi_deg2pix(degree, sys)
%deg2pix(degree, sys)
% Calculates pixel position from visual angle (degree) in a hemi-field (see
% picture). You can give it the screen 'width' or 'height' in a field of 
% sys or specify a window pointer field 'win' to let PTB estimate it.
% If they are provided as fields these are always used, irrespective 
% of whether they can be estimated from 'win' as well. Another field 'res'
% needs to include the screen resolution in pixel (e.g. [1680 1050]).
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
%   sys.dist     = 50; % cm
%   eccentricity = [10 0]; % degrees
%   pix = dfi_deg2pix(eccentricity, sys)
%
% ----------------------------------------
% adapted from Agoston Mihalic, 
% last updated, April 2015
%


    % enough information?
    if ~isfield(sys, 'width') || ~isfield(sys, 'height')
        if ~isfield(sys, 'win')
            error('I do not know the size of your monitor. Please provide it or a window pointer for estimation.')
        else
            disp('Estimating display size, because you did not provide height and width. This may be inaccurate.')
            [width, height]=Screen('DisplaySize', sys.win);
            width = width/10; % convert to cm 
            height= height/10;
        end
    end
    if isfield(sys, 'width'),  width  = sys.width;  end;
    if isfield(sys, 'height'), height = sys.height; end;
    
    % convert degrees of visual angle to pixels (along x and y)
    pixel = nan(size(degree));
    if exist('width', 'var')
        ppcmw      = sys.res(1)/width;
        pixel(:,1) = tand(degree(:,1)) * sys.dist * ppcmw;
    end
    % height
    if exist('height', 'var')
        ppcmh      = sys.res(2)/height;
        pixel(:,2) = tand(degree(:,2)) * sys.dist * ppcmh;
    end
    
end % end deg2pix

% eof


















