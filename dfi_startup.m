
% Calling dfi_startup only makes sense when in the right folder 
if ~exist('dfi', 'dir')
    warning('dfi_startup failed... go to dfi parent folder before calling this function.')
    return
end

% Call Psychtoolbox-3 specific startup function:
if exist('PsychStartup', 'file')
    fprintf('Initializing Psychtoolbox: Running PsychStartup...\n');
    PsychStartup; 
else
    warning('File "PsychStartup" not found!');
end

% % add fieldtrip folder
% fieldtripDir = fullfile('fieldtrip-20170929', 'fieldtrip-master');
% if exist(fieldtripDir, 'dir')
%     fprintf('Adding FieldTrip folder "%s" to Matlab search path ...\n', fieldtripDir);
%     addpath(fullfile(cd, fieldtripDir));
%     addpath(fullfile(cd, fieldtripDir, 'fileio'))
%     ft_defaults;
% else
%     warning(['Could not find FieldTrip folder "' , fieldtripDir, '"!'])
% end;

% add eeglab folder
eeglabDir = 'eeglab14_1_1b';
if exist(eeglabDir, 'dir')
    fprintf('Adding eeglab folder "%s" to Matlab search path ...\n', eeglabDir);
    %addpath(genpath(eeglabDir));
    addpath(fullfile(cd, eeglabDir));
else
    warning(['Could not find eeglab folder "' , eeglabDir, '"!'])
end

% add palamedes folder
palamedesDir = 'Palamedes';
if exist(palamedesDir, 'dir')
    fprintf('Adding Palamedes folder "%s" to Matlab search path ...\n', palamedesDir);
    addpath(fullfile(cd, palamedesDir));
else
    warning(['Could not find Palamedes folder "' , palamedesDir, '"!'])
end

% add SPM folder
spmDir = 'spm12';
if exist(spmDir , 'dir')
    fprintf('Adding SPM folder "%s" to Matlab search path ...\n', spmDir );
    addpath(fullfile(cd, spmDir));
else
    warning(['Could not find SPM folder "' , spmDir , '"!'])
end


% Make sure Matlab figures are saved with the background color as shown in
% Matlab figure gui
set(0, 'DefaultFigureInvertHardCopy', 'off');


% Re-seed random number generator using current time
rng('default')
rng('shuffle')

