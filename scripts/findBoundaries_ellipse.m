function [R1, R2, alpha_vec] = findBoundaries_ellipse(reference, ...
    curvesX, curvesY, curvesZ, T, N, B)
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

neighborhood = 300;

%Get dimensions
point_num = size(reference, 1);
trajectory_num = size(curvesX, 2);

R1 = zeros(point_num, 1);
R2 = zeros(point_num, 1);
alpha_vec = zeros(point_num, 1);
for ii = 1:point_num %For each point
    
    %Current Point
    point1 = reference(ii,:);
    
    %Get indices around point
    idx1 = ii-neighborhood;
    idx1 = max(idx1, 1);
    idx2 = ii+neighborhood;
    idx2 = min(idx2, point_num);  
    
    distances = zeros(trajectory_num, 1);
    max_point = zeros(trajectory_num, 3);
    for jj = 1:trajectory_num %For each trajectory     
        
        %Closest point to plane
        min_func = T(1,ii)*(curvesX(idx1:idx2,jj) - point1(1)) + ...
            T(2,ii)*(curvesY(idx1:idx2,jj) - point1(2)) + ...
            T(3,ii)*(curvesZ(idx1:idx2,jj) - point1(3));
        [~, ind2] = min(abs(min_func));
        point2 = [curvesX(idx1-1+ind2, jj), curvesY(idx1-1+ind2, jj), ...
            curvesZ(idx1-1+ind2, jj)];
        
        %Store temporary values
        max_point(jj, :) = point2;
        distances(jj) = norm(point1-point2);
    end
    
    for k=1:size(max_point,1)
        new_points(k,:) = [dot(N(:,ii), max_point(k,:)-point1), dot(B(:,ii), max_point(k,:)-point1)];
    end
    new_center = [0,0];
    [a,b,vect_a,vect_b,alpha] = estimateEllipse(new_points', new_center');
    
    R1(ii) = 1.01*a;
    R2(ii) = 1.01*b;
    alpha_vec(ii) = alpha;
    
end
end