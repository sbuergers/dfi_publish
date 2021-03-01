function stim = dfi_prepViStim(s, data)
% prepare blob/disk stimuli for fast screen presentation
%
% ----------------------------------------
% adapted from Agoston Mihalik, 
% last updated, April 2015
%
    disp    = s.disp;
    stim    = s.stim;
    ntrials = size(data, 1); 
    diam    = dfi_deg2pix(stim.vdiam, disp);
    dist    = dfi_deg2pix(stim.vloc, disp);
    
    % create stimuli with varying contrast per trials if requested
    numtrls = s.ntrls;
    if ~stim.vis.jittercntr && ~stim.vis.jitterpos, numtrls = 1; end
    
    % create contrast for each trial
    if isfield(stim.vis, 'contrbnds')
        contrbnds = stim.vis.contrbnds;
    else
        contrbnds = [0.1 0.9];
    end
%     Generate values from the uniform distribution on the interval [a, b].
%     r = a + (b-a).*rand(100,1);
    contrarr = contrbnds(1) + (contrbnds(2)-contrbnds(1)) .* rand(numtrls, 1);        
    if ~stim.vis.jittercntr
        contrarr = stim.vis.contr;
        numtrls  = 1;
    end
    
    % trial-loop
    for itrl = 1:numtrls
        
        % current trial's contrast
        s.stim.vis.contr = contrarr(itrl);
        
        % single trial stimuli
        if strcmp(stim.vis.type, 'circle')
            stim.vis.rect = CenterRectOnPoint([0 0 diam(1), diam(2)], disp.xcen-dist(1), disp.ycen+dist(2));
        elseif strcmp(stim.vis.type, 'blob')            
            stim.vis.image    = dfi_generate_gauss_blob(s);
            stim.vis.texture  = Screen(disp.win, 'MakeTexture', stim.vis.image);
            stim.vis.rect     = CenterRectOnPoint([0 0 diam(1) diam(2)], disp.xcen-dist(1), disp.ycen+dist(2));
        elseif strcmp(stim.vis.type, '2blobs')            
            stim.vis.image    = dfi_generate_gauss_blob(s);
            stim.vis.texture  = Screen(disp.win, 'MakeTexture', stim.vis.image);
            stim.vis.texture2 = Screen(disp.win, 'MakeTexture', stim.vis.image);
            stim.vis.rect     = CenterRectOnPoint([0 0 diam(1) diam(2)], disp.xcen-dist(1), disp.ycen+dist(2));
            stim.vis.rect2    = CenterRectOnPoint([0 0 diam(1) diam(2)], disp.xcen+dist(1), disp.ycen+dist(2));
        elseif strcmp(stim.vis.type, 'test')
            stim.vis.image    = repmat(255, disp.res(2), disp.res(1));
            stim.vis.texture  = Screen(disp.win, 'MakeTexture', stim.vis.image);
            stim.vis.rect     = CenterRectOnPoint([0 0 disp.res(1) disp.res(2)], disp.xcen, disp.ycen);
        else
            error('Invalid graphics parameter...')
        end
        
        % include positional jitter?
        if stim.vis.jitterpos
            if isfield(stim.vis, 'jitperc')
                maxjit = stim.vis.jitperc; % maximum percentage of jitter relative to the stimulus diameter
            else
                maxjit = 0.2;
            end
            possdir = [-1 1];
            actdir  = Sample(possdir);
            jitx    = rand * maxjit * diam(1); jitx = jitx*actdir; 
            jity    = rand * maxjit * diam(2); jity = jity*actdir;
            stim.vis.rect = CenterRectOnPoint([0 0 diam(1) diam(2)], disp.xcen-dist(1)+jitx, disp.ycen+dist(2)+jity);
        end
        
        % add trial-trial varying images to structure?
        if stim.vis.jittercntr
            stim.vis.imarr{itrl}   = stim.vis.image;
            stim.vis.textarr{itrl} = stim.vis.texture;
            stim.vis.rectarr{itrl} = stim.vis.rect;
        end
    end
        
    
end % end ve_prepViStim



%% old
% function stim = dfi_prepViStim(s, data)
% % prepare blob/disk stimuli for fast screen presentation
% %
% % ----------------------------------------
% % adapted from Agoston Mihalic, 
% % last updated, April 2015
% %
%     disp    = s.disp;
%     stim    = s.stim;
%     ntrials = size(data, 1); 
%     diam    = dfi_deg2pix(stim.vdiam, disp);
%     dist    = dfi_deg2pix(stim.vloc, disp);
%     
%     
%     if strcmp(stim.vis.type, 'circle')
%         stim.vis.rect = CenterRectOnPoint([0 0 diam(1), diam(2)], disp.xcen-dist(1), disp.ycen+dist(2));
%     elseif strcmp(stim.vis.type, 'blob')            
%         stim.vis.image   = dfi_generate_gauss_blob(s);
%         stim.vis.texture = Screen(disp.win, 'MakeTexture', stim.vis.image);
%         stim.vis.rect    = CenterRectOnPoint([0 0 diam(1) diam(2)], disp.xcen-dist(1), disp.ycen+dist(2));
%     elseif strcmp(stim.vis.type, '2blobs')            
%         stim.vis.image    = dfi_generate_gauss_blob(s);
%         stim.vis.texture  = Screen(disp.win, 'MakeTexture', stim.vis.image);
%         stim.vis.texture2 = Screen(disp.win, 'MakeTexture', stim.vis.image);
%         stim.vis.rect     = CenterRectOnPoint([0 0 diam(1) diam(2)], disp.xcen-dist(1), disp.ycen+dist(2));
%         stim.vis.rect2    = CenterRectOnPoint([0 0 diam(1) diam(2)], disp.xcen+dist(1), disp.ycen+dist(2));
%     elseif strcmp(stim.vis.type, 'test')
%         stim.vis.image   = repmat(255, disp.res(2), disp.res(1));
%         stim.vis.texture = Screen(disp.win, 'MakeTexture', stim.vis.image);
%         stim.vis.rect    = CenterRectOnPoint([0 0 disp.res(1) disp.res(2)], disp.xcen, disp.ycen);
%     else
%         error('Invalid graphics parameter...')
%     end
%     
%     
% end % end ve_prepViStim




% eof