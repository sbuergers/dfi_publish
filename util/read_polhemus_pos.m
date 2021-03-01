function [elec] = read_polhemus_pos(filename);

% Adapted from read_brainvision_pos (see blow) by Steffen Buergers,
% the difference between this one and the original is that in my .pos data
% I have an additional column at the very beginning that indexes the
% row/the recorded point. This is unnecessary, but obviously has to be
% taken into account with this function, otherwise there will be a an
% assignment mismatch further down the road.
% sbuergers@gmail.com

% READ_BRAINVISION_POS reads electrode positions measured with the Polhemus
% tracker in one of the F.C. Donders EEG labs. The polhemus software is actually 
% not from Brainvision.
%
% Use as:
%   [elec] = read_brainvision_pos(filename)
%
% This returns an electrode structure with
%   elec.label     cell-array with electrode labels (strings)
%   elec.pnt       position of each electrode

% Copyright (C) 2004, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

fid = fopen(filename, 'rt');
line = fgetl(fid);
Nchan = str2double(line);

if Nchan > 66,
    Nchan = 66;
    fprintf('\n\nReducing number of channels to 66 (64 + ref and ground).\nWe are not interested in pre-auricular points here.\n\n');
end

for i=1:Nchan
  line = fgetl(fid);
  [t, r] = strtok(line);
  t      = strtok(r);
  elec.label{i} = char(t);
  if ~isnan(sscanf(r, '%f'))
    tmp = sscanf(r, '%f')';
    elec.pnt(i,:) = tmp(2:4);
  else
    elec.pnt(i,:) = sscanf(r(regexp(r, '[^/A-Z]')), '%f')';
  end
end
elec.label = elec.label(:);

try
  % read the fiducials
  % find first line where the first token is not a number
  
  line = fgetl(fid);
  [t, r] = strtok(line);
  [t, r] = strtok(r);
  fiducial.label{1} = char(t);
  fiducial.pnt(1,:) = sscanf(r, '%f')';
  line = fgetl(fid);
  [t, r] = strtok(line);
  [t, r] = strtok(r);
  fiducial.label{2} = char(t);
  fiducial.pnt(2,:) = sscanf(r, '%f')';
  line = fgetl(fid);
  [t, r] = strtok(line);
  [t, r] = strtok(r);
  fiducial.label{3} = char(t);
  fiducial.pnt(3,:) = sscanf(r, '%f')';
  % add the fiducials to the electrode array
  elec.label = cat(1, elec.label(:), fiducial.label(:));
  elec.pnt   = cat(1, elec.pnt, fiducial.pnt);
catch
  % do nothing
  fprintf('\n\nInclusion of fiducials did not work...\n\n')
end

fclose(fid);
