function ab = dfi_conf_butterfly_dave(s)
% Create a confidence judgment butterfly response screen
%
% Input: s     = structure array containing all relevant variables for the
%                PTB screen and the stimuli in the sub-arrays s.disp and
%                s.stim. The optional field s.butterfly has to contain all
%                specifications used to create the butterfly confidence
%                buttons on the screen. If it is omitted default values
%                will be used (specified at the very top of the script).
% Output: ab   = Structure that contains all relevant drawing parameters 
%                for the butterfly response buttons.
%
% Details:
%       The function draws a butterfly response window on the screen, which
%       contains the set number of buttons (default is four) on the left
%       and right side of the screen around a central circle called
%       'circle_in' with the spread from the center defined by
%       'circle_out'. The response buttons are meant to be used to record
%       confidence judgments together with a binary response (e.g. left or
%       right would be quite natural). 
%
% ----------------------------------------
% adapted from David Meijer, 
% last updated, March 2016
%

% load default values
if ~isfield(s, 'butterfly')
    butterfly.circle_in.size  = 1.5; % inner circle size in visual angle
    butterfly.circle_out.size = 8;   % outer circle size 
    butterfly.button_angle    = 18;  % width of a button in visual angle
    butterfly.button_number   = 4;   % unilateral number of response buttons
else
    butterfly = s.butterfly;
end

% initialize variables for easier handling
disp       = s.disp; 
stim       = s.stim;
b_angle    = butterfly.button_angle;
b_num      = butterfly.button_number;
circle_in  = butterfly.circle_in;
circle_out = butterfly.circle_out;
ab.b_angle = b_angle;

% get screen distance
if ~isfield(disp, 'width') || ~isfield(disp, 'height')
    error('You need to provide the fields disp.width and disp.height for drawing the confidence response screen.');
end
if isfield(disp, 'width'),  width  = disp.width;  end;
if isfield(disp, 'height'), height = disp.height; end;

% DEBUG:
disp.res    = [1920 1080]
disp.width  =  48

    
% prepare inner circle (convert angle to rectangle of pixels)
circle_in.diam      = circle_in.size; 
pix_widths         = dfi_deg2pix([circle_in.diam, circle_in.diam], disp); 
circle_in.pix_wdth = pix_widths(1);         
circle_in.pix_hght = pix_widths(2); 
circle_in.rect     = CenterRectOnPoint([0 0 pix_widths(1) pix_widths(2)], disp.xcen, disp.ycen);
ab.inner_rect      = circle_in.rect;
clear pix_widths

% prepare outer circle (convert angle to rectangle of pixels)
circle_out.diam = circle_out.size;        
pix_widths     = dfi_deg2pix([circle_out.diam, circle_out.diam], disp);
ab.x_pix       = pix_widths(1);
ab.y_pix       = pix_widths(2); 
ab.outer_rect  = CenterRectOnPoint([0 0 pix_widths(1) pix_widths(2)], disp.xcen, disp.ycen);
clear pix_widths   

% angles of button lines to be drawn    
ab.part_angles    = [b_angle*2, b_angle, 0, ...
                    -b_angle,  -b_angle*2];
ab.part_angles_r  = 90  + ab.part_angles;                                                                        
ab.part_angles_l  = 270 - ab.part_angles;                                                                       
ab.full_angles_r  = ab.part_angles_r([1 5]);
ab.full_angles_l  = ab.part_angles_l([1 5]);


% get height and width of intersection point with the inner circle
ab.inner_x_angl = cosd(ab.part_angles) * circle_in.diam;
ab.inner_y_angl = sind(ab.part_angles) * circle_in.diam;
% transmute visual angles to pixels
pixels = dfi_deg2pix([ab.inner_x_angl', ab.inner_y_angl'], disp);
% and apply to left and right side of the display
ab.inner_x_pxls_r = disp.xcen + pixels(:,1)./2;
ab.inner_x_pxls_l = disp.xcen - pixels(:,1)./2;
ab.inner_y_pxls_r = disp.ycen + pixels(:,2)./2;
ab.inner_y_pxls_l = disp.ycen - pixels(:,2)./2;
clear pixels

% now do the same for the outer circle of the response buttons
ab.outer_x_angl = cosd(ab.part_angles)*circle_out.diam;
ab.outer_y_angl = sind(ab.part_angles)*circle_out.diam;
pixels = dfi_deg2pix([ab.outer_x_angl', ab.outer_y_angl'], disp);
ab.outer_x_pxls_r = disp.xcen + pixels(:,1)./2;
ab.outer_x_pxls_l = disp.xcen - pixels(:,1)./2;
ab.outer_y_pxls_r = disp.ycen + pixels(:,2)./2;
ab.outer_y_pxls_l = disp.ycen - pixels(:,2)./2;
clear pixels
  

% Prepare button lines, format data for the DrawLines Screen-function:
% first row is x coordinates, second row is y coordinates.
% two consecutive columns form a pair of (x,y) coordinates (start, end).
line_parts_r = [];
line_parts_l = [];
line_full_r  = [];
line_full_l  = [];
for i=1:5
    line_R_start = [ab.inner_x_pxls_r(i); ab.inner_y_pxls_r(i)];
    line_R_end   = [ab.outer_x_pxls_r(i); ab.outer_y_pxls_r(i)];
    line_parts_r = [line_parts_r line_R_start line_R_end];

    line_L_start = [ab.inner_x_pxls_l(i); ab.inner_y_pxls_l(i)];
    line_L_end   = [ab.outer_x_pxls_l(i); ab.outer_y_pxls_l(i)];
    line_parts_l = [line_parts_l line_L_start line_L_end];

    if i == 1 || i == 5
        line_full_r = [line_full_r line_R_start line_R_end];
        line_full_l = [line_full_l line_L_start line_L_end];
    end
end

% Draw on screen
% outer circle
Screen('FrameArc', w, grey, ab.outer_rect, ab.full_angles_r(2),(b_angle * b_num),3);
Screen('FrameArc', w, grey, ab.outer_rect, ab.full_angles_l(1),(b_angle * b_num),3);
% inner circle
Screen('FrameArc', w, grey, ab.inner_rect, ab.full_angles_r(2),(b_angle * b_num),3);
Screen('FrameArc', w, grey, ab.inner_rect, ab.full_angles_l(1),(b_angle * b_num),3);
% button lines
Screen('DrawLines', w, line_parts_r, 3, grey);
Screen('DrawLines', w, line_parts_l, 3, grey);
Screen(w, 'Flip');

% the lower and upper answering boxes
full_angl_bottom     = 180 + ab.part_angles([5 1]);   % relative to yAxis (clockwise)
full_angl_bottom_tmp =  90 - ab.part_angles([1 5]);   % relative to xAxis (anti-clockwise)
full_angl_top        =   0 + ab.part_angles([5 1]);   % relative to yAxis (clockwise)
full_angl_top_tmp    = -90 - ab.part_angles([1 5]);   % relative to xAxis (anti-clockwise)

inner_x_angl_bottom  = sind(full_angl_bottom_tmp)*circle_in.diam;
inner_x_pixls_bottom = disp.ycen + tand(inner_x_angl_bottom)*P.dist2Screen_nPixHeight;
ab.inner_x_angl_Top = sind(full_angl_top_tmp)*circle_in.diam;
ab.innerHeights_pixels_Top = disp.ycen + tand(ab.inner_x_angl_Top)*P.dist2Screen_nPixHeight;

ab.inner_y_angl_Bottom  = cosd(full_angl_bottom_tmp)*circle_in.diam;
ab.innerWidths_pixels_Bottom  = disp.xcen + tand(ab.inner_y_angl_Bottom)*P.dist2Screen_nPixWidth;
ab.inner_y_angl_Top  = cosd(full_angl_top_tmp)*circle_in.diam;
ab.innerWidths_pixels_Top  = disp.xcen + tand(ab.inner_y_angl_Top)*P.dist2Screen_nPixWidth;

inner_x_angl_bottom = sind(full_angl_bottom_tmp)*circle_out.diam;
ab.outerHeights_pixels_Bottom = disp.ycen + tand(inner_x_angl_bottom)*P.dist2Screen_nPixHeight;
ab.inner_x_angl_Top = sind(full_angl_top_tmp)*circle_out.diam;
ab.outerHeights_pixels_Top = disp.ycen + tand(ab.inner_x_angl_Top)*P.dist2Screen_nPixHeight;

ab.inner_y_angl_Bottom  = cosd(full_angl_bottom_tmp)*circle_out.diam;
ab.outerWidths_pixels_Bottom  = disp.xcen + tand(ab.inner_y_angl_Bottom)*P.dist2Screen_nPixWidth;
ab.inner_y_angl_Top  = cosd(full_angl_top_tmp)*circle_out.diam;
ab.outerWidths_pixels_Top  = disp.xcen + tand(ab.inner_y_angl_Top)*P.dist2Screen_nPixWidth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% confidence scale (double arrow lines) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ab.scale_angle = circle_out.diam*2;
ab.scale_Horpix = round(tand(ab.scale_angle)*P.dist2Screen_nPixWidth);          
ab.scale_highest_point = ab.inner_x_pxls_l(1);
ab.scale_lowest_point = ab.inner_x_pxls_l(5);
ab.scale_size_arrow = (ab.scale_lowest_point - ab.scale_highest_point)*(1/8);

ab.scale_X1 = ab.scale_Horpix - ab.scale_size_arrow;
ab.scale_X2 = ab.scale_Horpix;
ab.scale_X3 = ab.scale_Horpix + ab.scale_size_arrow;
ab.scale_X_coords = [ab.scale_X2, ab.scale_X2, ab.scale_X1, ab.scale_X2, ab.scale_X2, ab.scale_X3, ab.scale_X1, ab.scale_X2, ab.scale_X2, ab.scale_X3];
ab.scale_X_coords_R = disp.xcen + ab.scale_X_coords;
ab.scale_X_coords_L = disp.xcen - ab.scale_X_coords;

ab.scale_Y1 = ab.scale_lowest_point;
ab.scale_Y2 = ab.scale_lowest_point-ab.scale_size_arrow;
ab.scale_Y3 = ab.scale_highest_point+ab.scale_size_arrow;
ab.scale_Y4 = ab.scale_highest_point;
ab.scale_Y_coords = [ab.scale_Y4, ab.scale_Y1, ab.scale_Y2, ab.scale_Y1, ab.scale_Y1, ab.scale_Y2, ab.scale_Y3, ab.scale_Y4, ab.scale_Y4, ab.scale_Y3];



% eof






