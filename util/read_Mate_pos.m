function [sens, shape] = read_Mate_pos(filename)

% this function reads channels and fiducials positions from Polhemus file,
% which was created based on Mate's template. Fieltrip has couple other
% functions that read pos files (read_brainvision_pos & read_ctf_pos), but
% they are not compatibile with our Polhemus files. Output can be used to
% create layouts.

fid = fopen(filename, 'rt');
line = fgetl(fid);
Nchan = str2double(line);

if isnan(Nchan)
    Nchan = 69;
    fclose(fid);
    fid = fopen(filename, 'rt');
end

for i=1:Nchan
  line = fgetl(fid);
  C(i,:)= strsplit(line); 
end

try
    elec.pnt(1:64,1:3) = str2double(C(1:64,3:5)); % electrode positions
catch
    elec.pnt(1:64,1:3) = str2double(C(1:64,2:4)); % electrode positions
end

try
    for i=1:6
        line = fgetl(fid);
        Cf(i,:)= strsplit(line);
    end
    % calculate avg fiducials positions

    NA(1,1:3)=Cf(1,2:4);
    NA(2,1:3)=Cf(4,2:4);
    NA(3,1:3)=C(67,3:5);
    NA=str2double(NA);
    fiducial.pnt(1,:)=mean(NA,1);

    LPA(1,1:3)=Cf(2,2:4);
    LPA(2,1:3)=Cf(5,2:4);
    LPA(3,1:3)=C(68,3:5);
    LPA=str2double(LPA);
    fiducial.pnt(2,:)=mean(LPA,1);

    RPA(1,1:3)=Cf(3,2:4);
    RPA(2,1:3)=Cf(6,2:4);
    RPA(3,1:3)=C(69,3:5);
    RPA=str2double(RPA);
    fiducial.pnt(3,:)=mean(RPA,1);
catch
    fiducial.pnt(1,:) = str2double(C(67,2:4));
    fiducial.pnt(2,:) = str2double(C(68,2:4));
    fiducial.pnt(3,:) = str2double(C(69,2:4));
end

% add fiducials
elec.pnt   = cat(1, elec.pnt, fiducial.pnt);

%assign labels
elec.label={'Fp1'; 'Fp2'; 'F7'; 'F3'; 'Fz'; 'F4'; 'F8'; 'FC5'; 'FC1'; 'FC2'; 'FC6'; 'T7'; 'C3'; 'Cz'; 'C4'; 'T8';...
            'TP9'; 'CP5'; 'CP1'; 'CP2'; 'CP6'; 'TP10'; 'P7'; 'P3'; 'Pz'; 'P4'; 'P8'; 'PO9'; 'O1'; 'Oz'; 'O2'; 'PO10';...
            'AF7'; 'AF3'; 'AF4'; 'AF8'; 'F5'; 'F1'; 'F2'; 'F6'; 'FT9'; 'FT7'; 'FC3'; 'FC4'; 'FT8'; 'FT10'; 'C5'; 'C1';...
            'C2'; 'C6'; 'TP7'; 'CP3'; 'CPz'; 'CP4'; 'TP8'; 'P5'; 'P1'; 'P2'; 'P6'; 'PO7'; 'PO3'; 'POz'; 'PO4'; 'PO8';...
            'NAS'; 'LPA'; 'RPA'};
fiducial.label= {'NAS'; 'LPA'; 'RPA'};

fclose(fid);

shape.pnt=elec.pnt; % all positions for headshape
shape.fid.pnt=fiducial.pnt; % fiducials positions
shape.fid.label=fiducial.label; % fiducials labels
shape = ft_struct2double(shape);
sens = ft_datatype_sens(elec);

% eof
























