function [R1, R2, alpha_vec, vect_a_vec, vect_b_vec] = find_boundaries_ellipse(reference, ...
    curvesX, curvesY, curvesZ)
%Get dimensions
point_num = size(reference, 1);
trajectory_num = size(curvesX, 2);

R1 = zeros(point_num, 1);
R2 = zeros(point_num, 1);
alpha_vec = zeros(point_num, 1);
vect_a_vec = zeros(point_num, 3);
vect_b_vec = zeros(point_num, 3);
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
    
    vect_a = vect_a/norm(vect_a); vect_b = vect_b/norm(vect_b);
    project_coordinates = [size(points, 1), 2];
    for jj = 1:size(points,1)
       project_coordinates(jj, 1) = dot(points(jj,:) - means, vect_a);
       project_coordinates(jj, 2) = dot(points(jj,:) - means, vect_b);
    end
    
    %Check each magnitude of minor axis from smallest to largest
    %Finish when all points are enclosed by ellipse
    done_flag = false;
    for kk = 1:trajectory_num
        b = sorted_distances(kk);
        done_flag = check_ellipse(project_coordinates, a, b, alpha);
        if done_flag
            break;
        end
    end
    
    R1(ii) = a;
    R2(ii) = b;
    alpha_vec(ii) = alpha;
    vect_a_vec(ii,:) = vect_a/norm(vect_a);
    vect_b_vec(ii,:) = vect_b/norm(vect_b);
end
end

function flag = check_ellipse(points, a, b, alpha)
xp = points(:, 1); yp = points(:, 2);
val1 = (cos(alpha)*(xp) + sin(alpha)*(yp)).^2;
val2 = (sin(alpha)*(xp) - cos(alpha)*(yp)).^2;
if (val1/a^2)+(val2/b^2) <=1
    flag = true;
else
    flag = false;
end
end