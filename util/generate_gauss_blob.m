function I = generate_gauss_blob(setup, stim, dummy)
% generate_gauss_blob(setup, stim, dummy)
% 
% USAGE: 
%   I = generate_gauss_blob(setup, stim, dummy)
% DETAILS: 
%   Creates a Gaussian shaped two-dimensional blob of specified shape.
% INPUT: 
%   stim.vis.xstd ---     x-axis 
%   stim.vis.ystd --- and y-axis standard deviation
%   setup.ppd        ---
% OUTPUT:
%   I: . 
%
% ----------------------------------------
% adapted from Agoston Mihalic, 
% last updated, April 2015
%
 
% Calculate std in pixel
stdx = floor(stim.vis.xstd * setup.ppd); % ppd means pixels per degree
stdy = floor(stim.vis.ystd * setup.ppd);

% Define a meshgrid within 99.7% CI
[X, Y] = meshgrid(-3*stdx:3*stdx, -3*stdy:3*stdy);

% Calculate 2D Gaussian blob (with a simplified formula - easier, than mvnpdf)
I = exp(-(X.^2 / (2*stdx^2) + Y.^2 / (2*stdy^2))); 

% % create multivariate normal distribution 
% n  = 100000;
% t1 = normrnd(0, stdx, n, 1);
% t2 = normrnd(0, stdy, n, 1);
% covm= COV(t1,t2);
% I2= mvnpdf([X(:) Y(:)],0,covm);
% %     contour(X,Y,reshape(I2, sqrt(numel(X)), sqrt(numel(Y))));
% %     contour(X,Y,I);
% I = reshape(I2, sqrt(numel(X)), sqrt(numel(Y)));

% Rescale and contrast image
I = (I-min(I(:))) / (max(I(:))-min(I(:))); 
I = I * stim.vis.contr + (1-stim.vis.contr);
% I((x/3*sigma).^2+(y/3*sigma).^2>1) = 1-contrast; % would be out of 99.7% CI resampling

% Rescale for PTB or plot/write out 2D image
if dummy
%     figure;
%     colormap('jet');
%     mesh(X, Y, I)
    figure;
    colormap('gray');
    imshow(I);
%     imwrite(I, 'gauss.tiff')
else
    I = round(I * 255);
end

% eof