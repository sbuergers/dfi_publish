function D = ft2spm(cfg,ftData)
% Converts fieldtrip data structure to spm format

% List of condition labels,
condLabels  = cfg.condLabels;
fileName = cfg.fileName;

D = spm_eeg_ft2spm(ftData,fileName);

for c = 1:numel(condLabels)
    D = conditions(D,c,condLabels{c});
end

if isfield(cfg,'electrode_pos_file')
    % Assign actual EEG sensor positions
    % Loading chanel position file
    %elec = ft_read_sens(cfg.electrode_pos_file);
    elec = cfg.electrode_pos_file;
    elec = ft_convert_units(elec, 'mm');
    % Replacing the number labels in the electrode definition with
    % channel names
    % loadning channel map
    %elec.label(1:64) = cfg.channel_map(1:64,3);
    % Making sure the channel positions are identical to the
    % electrode positions, otherwise the alignment step does not work
    %elec.chanpos = elec.elecpos;
    
    fid.fid = struct();
    fid.fid.label = elec.label(end-2:end);
    fid.fid.pnt = elec.chanpos(end-2:end,:);
    fid.unit = elec.unit;
    fid.pnt = [];
    
    elec.chanpos(65:end,:) = [];
    elec.chantype(65:end)  = [];
    elec.chanunit(65:end)  = [];
    elec.elecpos(65:end,:) = [];
    elec.label(65:end)     = [];
    
    D = sensors(D,'EEG',elec);
    D = fiducials(D,fid);
    
else
    % Assign default EEG sensor positions
    S = struct();
    S.task = 'defaulteegsens';
    S.updatehistory = 0;
    S.D = D;
    D = spm_eeg_prep(S);
end

% Create 2D positions for EEG by projecting the 3D positions to 2D
S = struct();
S.task = 'project3D';
S.modality = 'EEG';
S.updatehistory = 0;
S.D = D;
D = spm_eeg_prep(S);

% at this stage units are 'unknown', replace with 'uV' for all channels
D = units(D, D.indchantype('EEG'), 'uV');

D.save

end

