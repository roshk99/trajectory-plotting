function [Router, xyz_distance_vec] = ...
    findBoundaries_circle(reference, curvesX, curvesY, curvesZ, ...
    neighborhood)
% -----------------------------------------------------------------------
% A function that finds the radius of a circular cross-section covering a
% set of trajectories using a plane approximation
%
% Inputs:
%   reference: the mean trajectory
%   curvesX, curvesY, curvesZ: the x, y, and z components of the individual
%                              trajectories
%   neighborhood: +/- the number of points to check to find the radius
%
% Output:
%   Router: the outer radius of the circular cross section for each "t"
%           value
%   xyz_distance_vec: the distance from the outermost point on the circle
%                     to the mean in components
%
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------

%Get dimensions
point_num = size(reference, 1);
trajectory_num = size(curvesX, 2);

%Calculate slopes
dx = gradient(reference(:,1)); dy = gradient(reference(:,2));
dz = gradient(reference(:,3));

Router = zeros(point_num, 1);
xyz_distance_vec = zeros(point_num, 3);
for ii = 1:point_num %For each point
    
    %Normalize gradient
    mag = norm([dx(ii) dy(ii) dz(ii)]);
    dx(ii) = dx(ii)/mag; dy(ii) = dy(ii)/mag; dz(ii) = dz(ii)/mag;
    
    %Get indices around point
    idx1 = ii-neighborhood;
    idx1 = max(idx1, 1);
    idx2 = ii+neighborhood;
    idx2 = min(idx2, point_num);
    
    %Current Point
    point1 = reference(ii,:);
    
    distances = zeros(trajectory_num, 1);
    xyz_distance = zeros(trajectory_num, 3);
    max_point = zeros(trajectory_num, 3);
    for jj = 1:trajectory_num %For each trajectory

        %Closest point to plane
        min_func = dx(ii)*(curvesX(idx1:idx2,jj) - point1(1)) + ...
            dy(ii)*(curvesY(idx1:idx2,jj) - point1(2)) + ...
            dz(ii)*(curvesZ(idx1:idx2,jj) - point1(3));
        [~, ind2] = min(abs(min_func));
        point2 = [curvesX(idx1-1+ind2, jj), curvesY(idx1-1+ind2, jj), ...
            curvesZ(idx1-1+ind2, jj)];        
        
        %Store temporary values
        max_point(jj, :) = point2;
        distances(jj) = norm(point1 - point2);
        xyz_distance(jj, :) = abs(point1 - point2);
    end
    
    %Get the point with the max distance
    [max_val, max_ind] = max(distances);
    
    %Store values
    Router(ii) = 1.01*max_val;
    xyz_distance_vec(ii, :) = xyz_distance(max_ind, :);
end
end