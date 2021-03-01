function I = dfi_generate_gauss_blob(s)
% dfi_generate_gauss_blob(s)
%
% s should contain
% s.disp.res  - screen resolution
% s.disp.win  - window pointer 
% 
% ----------------------------------------
% adapted from Agoston Mihalic, 
% last updated, April 2015
%


 
% Calculate std in pixel
stim = s.stim;
disp = s.disp;
pixel= dfi_deg2pix(s.stim.vdiam, disp); % pixel = [4 4]
stdx = pixel(1)/3;
stdy = pixel(2)/3;

% Define a meshgrid within 99.7% CI
[X, Y] = meshgrid(-3*stdx:3*stdx, -3*stdy:3*stdy);

% Calculate 2D Gaussian blob with a simplified formula 
I = exp(-(X.^2 / (2*stdx^2) + Y.^2 / (2*stdy^2))); 
clear X Y

% Rescale and contrast image
I = (I-min(I(:))) / (max(I(:))-min(I(:))); 
if ~s.disp.incolor
    I   = I * (255 * stim.vis.bgcontr * stim.vis.contr);
    I   = I + 255 * stim.vis.bgcontr;
    I   = round(I);
else
    I   = I * stim.vis.contr;
    I   = round(I * 255);
    R   = I;
    G   = repmat(255, size(R));
    G   = G-R;
    B   = zeros(size(I));
    RGB = [R G B];
    RGBplane = reshape(RGB, [size(I) 3]);
    clear RGB R G B
    I   = RGBplane;
    clear RGBplane
end;



% eof