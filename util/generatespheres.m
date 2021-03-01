function spheres = generatespheres(stim, criteria, isplot)
%Generate and optionally plot two dimensional, filled spheres, whose
% visibility behaves as a Gaussian 2-dimensional function, with mean [0,0]
% and standard deviation given in stim.visual.xstd, stim.visual.ystd. Visibility
% is highest for the functions mode.
%
% ----------------------------------------
% adapted from Agoston Mihalik, 
% last updated, April 2015
%
      
% Default randomization for 2D Gaussian
mu = [0 0];
sigma = [stim.visual.xstd^2 0; 0 stim.visual.ystd^2];
spheres = mvnrnd(mu, sigma, stim.visual.num);

% Resample spheres
idout = (spheres(:,1) ./ (2*stim.visual.xstd)).^2 + (spheres(:,2) ./ (2*stim.visual.ystd)).^2 > 1; % out of the 2D gaussian's 95% confidence interval
if strcmp(criteria, '2-fold') % 2-fold ctireria
    while std(spheres(:,1)) > stim.visual.xstd*1.1 || std(spheres(:,1)) < stim.visual.xstd*0.9 || sum(idout) % or std is not close to the expected value
        spheres = mvnrnd(mu, sigma, stim.visual.num); % resample all the spheres
        
        % Resample only spheres that are out of the 2D gaussian's 95% confidence interval
        idout = (spheres(:,1) ./ (2*stim.visual.xstd)).^2 + (spheres(:,2) ./ (2*stim.visual.ystd)).^2 > 1;
        while sum(idout)
            spheres(idout,:) = mvnrnd(mu, sigma, sum(idout));
            idout = (spheres(:,1) ./ (2*stim.visual.xstd)).^2 + (spheres(:,2) ./ (2*stim.visual.ystd)).^2 > 1;
        end
    end
elseif strcmp(criteria, '1-fold') % 1-fold criteria
    while sum(idout)
        spheres(idout,:) = mvnrnd(mu, sigma, sum(idout));
        idout = (spheres(:,1) ./ (2*stim.visual.xstd)).^2 + (spheres(:,2) ./ (2*stim.visual.ystd)).^2 > 1;
    end
end

% Plot if needed
if isplot
    figure;
    N = 50;
    theta = 0:1/N:2*pi;
    x = cos(theta) * 2 * stim.visual.xstd; 
    y = sin(theta) * 2 * stim.visual.ystd;
    plot(x, y, 'r');
    hold on
    plot(spheres(:,1), spheres(:,2), 'b.');
    axis equal;
end
% eof