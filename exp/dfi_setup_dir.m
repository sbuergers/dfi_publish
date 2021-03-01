function mypath = dfi_setup_dir(s)
% mypath = dfi_setup_dir(subj, flag)
% 
% DESCRIPTION:
%   setup_dir creates a structure 'mypath' that contains the complete
%   pathnames of folders used for saving or analyzing subject data, e.g.:
%
%   mypath.root:    home folder of this function
%   mypath.misc:    root/miscellaneous
%   mypath.data:    root/data/<subj_id>/pilot/<date>
%   or              root/data/<subj_id>/<date>
%
%   Depending on the data in 's.subj' a folder for pilots or experimental
%   sessions is added including an id number. I have excluded 'flag' for
%   now.
%
% ----------------------------------------
% adapted from Agoston Mihalik, 
% last updated, April 2015
%

subj = s.subj;

% Root folder
pathname    = fileparts(which('dfi_setup_dir.m')); 
mypath.root = fullfile(pathname, filesep);

% Experiment folder
mypath.data = fullfile(mypath.root, 'data', s.paradigm);

% Pilot and session data folder
if subj.pilot
    mypath.data = fullfile(mypath.data, subj.id, 'pilot', s.cond);
else
    mypath.data = fullfile(mypath.data, subj.id, s.cond);
end
mypath.data = fullfile(mypath.data, datestr(now, 'yyyy_mm_dd'));

% Create if not already existent
if ~isdir(mypath.data)
    mkdir(mypath.data) 
end

% What we did
disp('Adding folders... ');
disp(mypath);

% eof
















