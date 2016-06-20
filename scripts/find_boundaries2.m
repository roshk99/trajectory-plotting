function [R1, xyz_distance_vec, vect_a_vec, vect_b_vec] = ...
    find_boundaries2(reference, curvesX, curvesY, curvesZ)

%Get dimensions
point_num = size(reference, 1);
trajectory_num = size(curvesX, 2);

R1 = zeros(point_num, 1);
vect_a_vec = zeros(point_num, 3);
vect_b_vec = zeros(point_num, 3);
xyz_distance_vec = zeros(point_num, 3);
for ii = 1:point_num %For each point
    
    %Calculate the mean point
    means = reference(ii,:);
    
    %Get the current set of points and calculate distances from the mean
    points = [curvesX(ii, :)' curvesY(ii,:)' curvesZ(ii,:)'];
    distances = zeros(trajectory_num, 1);
    for jj = 1:trajectory_num
        distances(jj) = norm(points(jj,:) - means);
    end
    
    %Sort the distances
    [sorted_distances, indx] = sort(distances, 'ascend');
    
    %Major axis radius is maximum distance
    a = sorted_distances(end);
    maj_dist_point = points(indx(end), :);
    
    %Get the major axis vector and angle of rotation
    vect_a = maj_dist_point - means;
    alpha = atan2(vect_a(2), vect_a(1));
    
    vect2 = vect_a;
    for jj = 1:size(points, 1)
        vect2 = cross(vect2, points(jj,:)-means);
    end
    n = vect2/norm(vect2);
    vect_b = cross(vect_a, n);  
        
    R1(ii) = a;
    xyz_distance_vec(ii,:) = abs(vect_a);
    vect_a_vec(ii,:) = vect_a/norm(vect_a);
    vect_b_vec(ii,:) = vect_b/norm(vect_b);
end
end