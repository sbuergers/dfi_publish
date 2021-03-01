function I = dfi_generate_circle(s)
% dfi_generate_circlce(s)
%
% s should contain
% s.disp.res  - screen resolution
% s.disp.win  - window pointer 
% s.stim.vloc - degrees from center [horizontal, vertical]
% 
%
 
% Define the image size (in pixels)
diam = dfi_deg2pix(s.stim.vdiam, s.disp);
rows = diam(1);
cols = diam(2);
radius = 20;
center = [25 25];  % In [X,Y] coordinates
% Make the circle
[xMat,yMat] = meshgrid(1:cols,1:rows);
distFromCenter = sqrt((xMat-center(1)).^2 + (yMat-center(2)).^2);
circleMat = distFromCenter<=radius;
figure, imshow(circleMat)

pixel= deg2pix(s.stim.vloc, s.disp);
xpos = pixel(1);
ypos = pixel(2);

% Define a meshgrid containing the circle
[X, Y] = meshgrid(1:xpos, 1:ypos);

% 




% eof