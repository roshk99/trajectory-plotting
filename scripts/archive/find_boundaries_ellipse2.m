function [R1, R2, alpha_vec] = find_boundaries_ellipse(reference, ...
    curvesX, curvesY, curvesZ, T, N)
% -----------------------------------------------------------------------
% A function that finds the two radii and angle of orientation for ellipses
% as cross sections
%
% Inputs:
%   reference: the mean trajectory
%   curvesX, curvesY, curvesZ: the x, y, and z components of the individual
%                              trajectories
%   T, B: part of the TNB reference frame for each data point
%
% Output:
%   R1, R2: the major and minor radii of the ellipse
%   alpha: the angle between the major axis and the N-axis
%
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------

%Get dimensions
point_num = size(reference, 1);
trajectory_num = size(curvesX, 2);

R1 = zeros(point_num, 1);
R2 = zeros(point_num, 1);
alpha_vec = zeros(point_num, 1);
for ii = 1:point_num %For each point
    
    %Current Point
    point1 = reference(ii,:);
    
    distances = zeros(trajectory_num, 1);
    max_point = zeros(trajectory_num, 3);
    for jj = 1:trajectory_num %For each trajectory
        
        %Closest point to plane
        min_func = T(1,ii)*(curvesX(:,jj) - point1(1)) + ...
            T(2,ii)*(curvesY(:,jj) - point1(2)) + ...
            T(3,ii)*(curvesZ(:,jj) - point1(3));
        [~, ind2] = min(abs(min_func));
        point2 = [curvesX(ind2, jj), curvesY(ind2, jj), curvesZ(ind2, jj)];
        
        %Store temporary values
        max_point(jj, :) = point2;
        distances(jj) = norm(point1-point2);
    end
    
    [sorted_distances, indx] = sort(distances, 'ascend');
    
    %Get the point with the max distance
    a = sorted_distances(end);
    b = sorted_distances(end-1);
    
    maj_dist_point = max_point(indx(end), :);
    
    %Get the major axis vector
    vect_a = maj_dist_point - point1;
    vect_a = vect_a/norm(vect_a);
    vect_b = cross(vect_a, T(:,ii));
    vect_b = vect_b/norm(vect_b);
    
    %Angle between the major axis and N-axis
    alpha = acos(dot(vect_a, N(:,ii)));
    
    R1(ii) = a;
    R2(ii) = b;
    alpha_vec(ii) = alpha;
    
end
end