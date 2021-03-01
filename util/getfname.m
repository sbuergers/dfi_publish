function fname = getfname(folder, strspec)
% fname = getfname(folder, spec)
%  getfname  lists all the file/folder names to a cell array that match the 
% string specifier. You can use wildcards (*).
%  Input:
%       folder  = folder name in which we want to list files/subfolders
%       strspec = string specifier/pattern for file/folder names
%  Output:
%       fname = cell array of file/folder names
%
% ----------------------------------------
% adapted from Agoston Mikhalic, 
% last updated, April 2015
%

% List the files in the folder based on the string specifier
files = dir(fullfile(folder, strspec));

% Put the file/folder names into a cell
fname = cell(size(files));
for i=1:length(files)
    fname{i} = files(i).name;
end