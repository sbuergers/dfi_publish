function [distances, minid] = dist_3d(locvect, point)
% takes a matrix of the format Nx3, where each column is a coordinate, and
% a point with 1x3. Calculates the euclidian distance between each point of
% the matrix and the point vector.
    x = locvect(:,1);
    y = locvect(:,2);
    z = locvect(:,3);
    xp = point(1); 
    yp = point(2);
    zp = point(3);
    distances = sqrt((x-xp).^2 + (y-yp).^2 + (z-zp).^2);
    [~, minid] = min(distances);
    
    
    
return