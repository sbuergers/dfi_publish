function [ d ] = smoothlyjointssegments( spec, dvect )
%d = SMOOTHLYJOINTSSEGMENTS(spec, dvect)
% takes two time series (in particular MEG or EEG data) segments (dvect contains data segments) 
% and concatenates them smoothly. The way the segments are smoothed
% around the "cut" is specified in spec.
%
% spec.Fs    = sampling frequency (default = 1000)
% spec.range = data points to be smoothed on each side (default = 50)
% spec.alpha = smoothing parameter of exponential moving average filter,
%              the average becomes smoother with smaller values for alpha (default = 0.4)
% spec.plot  = plot the original data over the adjusted data (default = 1)
% spec.plots = plot the whole range or only the modified segment (default =
%              0, i.e. only modified segment)
% spec.chan  = Channel to be plotted (default is 31)
% spec.numi  = Max number of images to be plotted (default = 5)
% spec.savei = save images or not (default = 0)
% spec.diri  = directory where the images should be saveid (default = current directory)
% spec.fn    = filename 
%
% OUTPUT: Gives a continous ts data segment after

if ~isfield(spec, 'Fs'), Fs = 1000;      else Fs    = spec.Fs;    end;
if ~isfield(spec, 'range'), range = 50;  else range = spec.range; end;
if ~isfield(spec, 'alpha'), alpha = 0.4; else alpha = spec.alpha; end;
if ~isfield(spec, 'mkplot'), mkplot = 1; else mkplot= spec.plot;  end;
if ~isfield(spec, 'plots'), plots = 0;   else plots = spec.plots; end;
if ~isfield(spec, 'chan'), ch = 31;      else ch    = spec.chan;  end;
if ~isfield(spec, 'numi'), numi = 5;     else numi  = spec.numi;  end;
if ~isfield(spec, 'savei'), savei = 0;   else savei = spec.savei; end;
if ~isfield(spec, 'diri'), diri = cd;    else diri  = spec.diri;  end;
if ~isfield(spec, 'fn'),fn='smoothplot'; else fn    = spec.fn;    end;

if numel(dvect) < 2, 
    warning('There is only one data segment. No changes are made.'); 
    d = dvect;
    return
end;

plotcount = 0;
while numel(dvect) > 1
    % update dvect until only one segment is left.
    d1 = dvect{1}; d2 = dvect{2};
    dvect(2) = [];
    d  = [d1, d2];
    % interpolate by filterting around the "cut"
    dinterpol   = d(:, length(d1) - range : length(d1) + range);
    movavg      = filter(alpha, [1, alpha-1], dinterpol, [], 2);
    % get the "residuals"/difference values for the time
    % period we want to smooth
    res = dinterpol - movavg;
    % apply smoothing and weight according to proximity to "cut"
    mu = range/2 + 1; % taken as the cut-point
    sd = range/5;
    pdv = repmat(normpdf(1:range + 1, mu, sd), 64, 1);
    % scale probabilities (at the cut w = 1)
    w = pdv*(1/pdv(1,mu));
    dsmooth = d(:, length(d1)-range/2:length(d1)+range/2) - res(:, mu:mu+range).*w;
    % plot 
    if mkplot && plotcount < numi 
        fig = figure('Color', [1 1 1]);
        hold on
        if plots
             h1 = plot(1:length(movavg(ch,:)), movavg(ch,:), 'b', ...
                  mu:(mu+range), dsmooth(ch,:), 'g', ...
                  1:length(movavg(ch,:)), dinterpol(ch,:), 'r', 'LineWidth',2);
        else
             h1 = plot(mu:(mu+range), movavg(ch,mu:(mu+range)), 'b', ...
                  mu:(mu+range), dsmooth(ch,:), 'g', ...
                  mu:(mu+range), dinterpol(ch,mu:(mu+range)), 'r', 'LineWidth',2);
        end
        legend(h1, 'movavg', 'weighted movavg', 'data')
        title(['Smoothed edge of two ts segments (channel ',num2str(ch),')'])
        set(findall(gcf,'-property','FontSize'),'FontSize',12)
        if savei
            cd(diri)
            saveas(fig, [fn, '_', int2str(numel(dvect))], 'png')
        end
    end
    plotcount = plotcount + 1;
    % now replace the "cut" segment with the smoothed one
    d(:, length(d1)-range/2:length(d1)+range/2) = dsmooth;
    dvect(1) = {d};    
end % eof









