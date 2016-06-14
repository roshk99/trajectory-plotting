function [Router, xyz_distance_vec] = find_boundaries(reference, ...
    curvesX, curvesY, curvesZ)
    %Reference is the trajectory from which you are drawing the
    %perpendicular
    %curvesX is point_num by trajectory_num
    
    %Get dimensions
    point_num = size(reference, 1);
    trajectory_num = size(curvesX, 2);
    
    %Calculate slopes
    dx = gradient(reference(:,1)); dy = gradient(reference(:,2)); 
    dz = gradient(reference(:,3));
    
    Router = zeros(point_num, 1);
    xyz_distance_vec = zeros(point_num, 3);
    for ii = 1:point_num %For each point
                    
        %Current Point
        point1 = reference(ii,:);
        
        distances = zeros(trajectory_num, 1);
        xyz_distance = zeros(trajectory_num, 3);
        for jj = 1:trajectory_num %For each trajectory
            
            %Closest point to plane
            min_func = dx(ii)*(curvesX(:,jj) - point1(1)) + ...
                       dy(ii)*(curvesY(:,jj) - point1(2)) + ...
                       dz(ii)*(curvesZ(:,jj) - point1(3));
            [~, ind2] = min(abs(min_func));
            point2 = [curvesX(ind2, jj), curvesY(ind2, jj), curvesZ(ind2, jj)];
            
            distances(jj) = sqrt(sum((point1 - point2).^2));
            xyz_distance(jj, :) = abs(point1 - point2);
        end
        [max_val, max_ind] = max(distances);
        Router(ii) = max_val;
        xyz_distance_vec(ii, :) = xyz_distance(max_ind, :);
    end
end