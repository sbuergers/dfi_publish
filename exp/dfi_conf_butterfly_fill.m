function button_status = dfi_conf_butterfly_fill( s, ab )
% Create a confidence judgment butterfly response screen button is filled
%
% Input: s     = structure array containing all relevant variables for the
%                PTB screen and the stimuli in the sub-arrays s.disp and
%                s.stim. The optional field s.butterfly has to contain all
%                specifications used to create the butterfly confidence
%                buttons on the screen. If it is omitted default values
%                will be used (specified at the very top of the script).
% Input: ab    = Structure that contains all relevant drawing parameters 
%                for the butterfly response buttons.
%
% Output: button_status = Which button (if any) was pressed?
%
% Details:
%       see dfi_conf_butterfly_draw for details.
%
% ----------------------------------------
% adapted from David Meijer, 
% last updated, March 2016
%

disp        = s.disp; 
stimColor   = disp.grey + disp.grey * s.stim.vis.contr;
resp_angle  = ab.resp_angle;
part_angles = -ab.part_angles;
r_angles    = ab.part_angles_r;
l_angles    = ab.part_angles_l;

% where is the button?
if ab.x < disp.xcen, left = 1; else left = 0; end;
if ab.y < disp.ycen, top  = 1; else top  = 0; end;

% which button exactly?
bnum = 0;
for i = 2:ab.b_num+1
    if resp_angle < part_angles(i) && bnum == 0;
        bnum = i - 1;
    end
end

% fill button
if bnum ~= 0
    if left
        button_status = -bnum;
        Screen('FillArc', disp.win, stimColor, ab.outer_rect, l_angles(bnum), ab.b_angle)
        Screen('FillArc', disp.win, disp.grey, ab.inner_rect, 180, 180)
    else
        button_status = +bnum;
        Screen('FillArc', disp.win, stimColor, ab.outer_rect, -r_angles(bnum + 1)+90, ab.b_angle);
        %Screen('FillArc', disp.win, stimColor, ab.outer_rect, r_angles(bnum + 1), ab.b_angle);
        Screen('FillArc', disp.win, disp.grey, ab.inner_rect, 0, 180)
    end
    % inner circle
    Screen('FrameArc', disp.win, stimColor, ab.inner_rect, ab.full_angles_r(2),(ab.b_angle * ab.b_num),3);
    Screen('FrameArc', disp.win, stimColor, ab.inner_rect, ab.full_angles_l(1),(ab.b_angle * ab.b_num),3);
else
    button_status = bnum;
end


return

% eof




























