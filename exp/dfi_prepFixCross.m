function fix = dfi_prepFixCross(s)
% prepare red & gray fix-cross stimuli for fast screen presentation
%
% ----------------------------------------
% adapted from Agoston Mihalic, 
% last updated, April 2015
%

    disp = s.disp;
    
    [fixWidth, lineWidth] = deal(0.7, 0.06); 
    if s.disp.incolor
        [image, rect]     = makefixationcross(fixWidth, lineWidth, disp.ppd, disp.green, [255 0 0]); % red fixation cross
    else
        if isfield(s.stim.vis, 'visib')
           [image, rect]  = makefixationcross(fixWidth, lineWidth, disp.ppd, disp.grey, ...
                                              disp.white * s.stim.vis.visib); % grey fixation cross
        else
           [image, rect]  = makefixationcross(fixWidth, lineWidth, disp.ppd, disp.grey, ...
                                              disp.grey + disp.grey * s.stim.vis.contr); % grey fixation cross
        end
    end
    fix.texture           = Screen(disp.win, 'MakeTexture', image);
    fix.rect              = CenterRectOnPoint(rect, disp.xcen, disp.ycen);
end % end dfi_prepFixCross


