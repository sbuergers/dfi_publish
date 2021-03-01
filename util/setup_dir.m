function mypath = setup_dir(subj, flag)
% mypath = setup_dir(subj, flag)
% 
% DESCRIPTION:
%   setup_dir creates a structure 'mypath' that contains the complete
%   pathnames of folders used for saving or analyzing subject data, e.g.:
%
%   mypath.root:    home folder of this function
%   mypath.misc:    root/miscellaneous
%   mypath.data:    root/data/<sessiontype_id>/<subj_id>/date
%   e.g.            root/data/pilot_2         /149      /2015_04_22
%
%   Depending on the data in 'subj' a folder for pilots or experimental
%   sessions is added including an id number. 'flag' is optional and used
%   for special additional folder settings for 'lab' and 'anal'.
%
%   'lab': ... to do ...
%
%   'anal': ... to do ...
%
%
% ----------------------------------------
% adapted from Agoston Mikhalic, 
% last updated, April 2015
%

if ~exist('flag', 'var')
   flag = ''; 
end

% Root folder
pathname    = fileparts(which('setup_dir.m')); 
mypath.root = fullfile(pathname, filesep);

% Add miscellaneous script folder to the path
mypath.misc = fullfile(mypath.root, 'miscellaneous');
if ~isdir(mypath.misc), mkdir(mypath.misc); end; 
addpath(mypath.misc); 

% Experiment folder
mypath.data = fullfile(mypath.root, 'data');

% Pilot and session data folder
folder = {'data'};
for i=1:length(folder)
    if isfield(subj, 'pilot')
        mypath.(folder{i}) = fullfile(mypath.(folder{i}), sprintf('pilot_%d', subj.pilot), subj.id);
    else
        mypath.(folder{i}) = fullfile(mypath.(folder{i}), subj.id);
    end
    if isfield(subj, 'session')
        session = getfname(mypath.data, '201*');
        if subj.session <= length(session) %&& strcmp(flag, 'interactive') % session already exists
            if ~strcmp(flag, 'anal') && ~strcmp(session{subj.session}, datestr(now, 'yyyy_mm_dd')) % session date is different!!
                error('The given session folder already exists with another date.');
            end
            mypath.(folder{i}) = fullfile(mypath.(folder{i}), session{subj.session});
        elseif subj.session > length(session) && strcmp(flag, 'anal')
            error('Session folder does not exists. Check the data files you want to analyse.');
        else
            mypath.(folder{i}) = fullfile(mypath.(folder{i}), datestr(now, 'yyyy_mm_dd'));
        end
    end
    if ~isdir(mypath.(folder{i}))
        mkdir(mypath.(folder{i})) 
    end
end

% Analysis folder if needed
if strcmp(flag, 'anal')
    plot2svgdir = fullfile(mypath.misc, 'plot2svg_20120915');
    if isdir(plot2svgdir)
        addpath(plot2svgdir) % to enable svg export
    end
    mypath.anal = fullfile(mypath.data, 'analysis');
    if ~isdir(mypath.anal)
        mkdir(mypath.anal)
    end
end

% Additional dropbox folders if needed
if  strcmp(flag, 'lab')
    mypath.dropbox = fullfile(mypath.root, '..', 'Dropbox', 'MATLAB');
    mypath.copydata = fullfile(mypath.dropbox, 'data', subj.id);
    if ~isdir(mypath.copydata)
        mkdir(mypath.copydata)
    end
end

% what we did
disp('Adding folders... ');
disp(mypath);

% eof
















