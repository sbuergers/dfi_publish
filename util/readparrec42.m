%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:  Reads Philips PAR/REC V4.2 files                        %
% Author:       Ramin S. Sahebjavaher                                   %
% Date:         June 2011, modified October 2012                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [magn, phas, info] = readparrec42( filename, display )
basefilename = cd; cd(filename); localfilename = ls;
fid1 = fopen(localfilename(3,:)); fid2 = fopen(localfilename(4,:),'r');
xPAR = textscan(fid1,'%s', 'delimiter', '\n'); xPAR = char(xPAR{1});
xREC = fread(fid2,'uint16');
fclose(fid1); fclose(fid2);
cd(basefilename);

% Calculate matrix size
N_slc=0; N_dyn=0; N_card=0; diffTE=0; PEV1=0; PEV2=0; PEV3=0;
if isempty( sscanf(xPAR(size(xPAR,1)-2,:),'%f') )
    NrLineToRead = size(xPAR,1)-3;
else
    NrLineToRead = size(xPAR,1)-2;
end
for cnt1=101:NrLineToRead
    param = sscanf(xPAR(cnt1,:),'%f');
    if (param(1) > N_slc)
        N_slc = N_slc + 1;
        rec_matrix  = [param(10) param(11) N_slc]; %x,y,z
    end
    if (param(3) > N_dyn)
        N_dyn = N_dyn + 1;
    end
    if (param(4) > N_card)
        N_card = N_card + 1;
    end
end

% Read data
cnt2=0; magn = zeros([rec_matrix N_dyn]); phas = zeros([rec_matrix N_dyn]); off_center = zeros(N_slc,3);
for cnt1=101:NrLineToRead
    param = sscanf(xPAR(cnt1,:),'%f');
    RI          = param(12);
    RS          = param(13);
    SS          = param(14);
    
    P_scan      = param(9);
    gap         = param(24);
    rez         = [param(29) param(30) param(23)]; %x,y,z
    Orientation = param(26);
    TE          = param(31);
    N_average   = param(35);
    
    if (param(5) == 0) % mag
        PV   = reshape(xREC((rec_matrix(1)*rec_matrix(2)*cnt2+1):(rec_matrix(1)*rec_matrix(2)*(cnt2+1))),rec_matrix(1:2));
        tmp  = (PV'.*RS+RI)./(RS.*SS); %x,y,x,dyn
        magn(:,:,param(1),param(3))  = tmp;
        if (display) 
            imagesc(magn(:,:,param(1),param(3)));
        end
    elseif (param(5) == 3) % phase
        PV   = reshape(xREC((rec_matrix(1)*rec_matrix(2)*cnt2+1):(rec_matrix(1)*rec_matrix(2)*(cnt2+1))),rec_matrix(1:2));
        tmp  = (PV'.*RS+RI)./(RS.*SS); %x,y,x,dyn
        phas(:,:,param(1),param(3)) = tmp;
        if (display) 
            imagesc(tmp,[-pi,pi]);
        end
    end
    if (display) 
        %pause(display)
        colormap(gray)
        drawnow
    end
    off_center(param(1),:) = [param(20) param(21) param(22)];  %based on slice
    cnt2 = cnt2+1;
end
  
info = {'Matrix (pix)', rec_matrix, 'Res_XYZ (mm)', rez , 'N_slices', N_slc, 'N_dynamics', N_dyn, 'N_card', N_card, ...
        'Orientation', Orientation, 'offcenter_XYZ (mm)',off_center, 'TE (ms)', TE,'N_averages', N_average};