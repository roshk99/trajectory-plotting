function trajectories = createTrajectory_ellipse(values, randomflag, points)
% -----------------------------------------------------------------------
% A function that takes smoothed data and generates new trajectories based
% on a canal surface approach
%
% Inputs:
%   values: contains parameters of the encoded demonstrations
%   randomflag: boolean depending on whether you want random start
%               points for trajectories to be generated
%   points: if randomflag is true, points is an integer with the number of
%           points to generate. if randomflag is false, points is a matrix
%           of dimension point_num x 3 that has the initial points for the
%           new trajectories
%
% Output:
%   trajectories: cell with each element point_numx3 vector of generated
%                 trajectories
%
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------

%THIS CODE IS NOT DONE AND STILL HAS SOME BUGS

%Get the initial mean point
mean_vec = [values(1).xmean, values(1).ymean, values(1).zmean];

%The number of trajectories
data_size = size(values(1).N, 1);

%Manually setting randomflag to true (haven't implemented the else portion)
if randomflag == false
    randomflag = true;
    points = 1;
end

%if randomly generated initial points are desired
if randomflag
    initPoints = zeros(points, 3);
    all_RR_vecs = zeros(points, 2);
    
    for ii = 1:points
        done_flag = false;
        while done_flag == false
            %Generate random point in R-R space and check if in ellipse
            RR_vec = [values(1).R1(1)*(-1+2*rand(1)) values(1).R2(1)*(-1+2*rand(1))];
            val = (RR_vec(1)/values(1).R1(1))^2 + (RR_vec(2)/values(1).R2(1))^2;
            if val < 1
                done_flag = true;
            end
        end
        %Rotate by alpha to get point in N-B space
        NB_vec = ([cos(values(1).alpha(1)) -sin(values(1).alpha(1)); ...
            sin(values(1).alpha(1)) cos(values(1).alpha(1))]*RR_vec')';
        all_RR_vecs(ii, :) = RR_vec;
        
        %Convert to 3D space
        initPoints(ii,:) = NB_vec(1)*values(1).N(1,:) + ...
            NB_vec(2)*values(1).B(1,:) + mean_vec(1,:);
    end
else
    %Add code here
end

%Compute the trajectory based on the output from the canal and initial
%points
trajectories = {};
for ii = 1:length(initPoints(:,1))
    trajectories{ii} = compute_trajectory(values, data_size, ...
        initPoints(ii,:), all_RR_vecs(ii,:));
end
end

function newTraj = compute_trajectory(values, data_size, initPoint, RR_vec)

%Gets the mean vector
mean_vec = [values(1).xmean, values(1).ymean, values(1).zmean];

%Find the ratio
dist1 = norm(RR_vec);
trial_distances = linspace(0, values(1).R1(1), 1000)';
trial_points = repmat(RR_vec/dist1, length(trial_distances), 1) ...
    .*repmat(trial_distances, 1, 2);
vals = (trial_points(:,1)./values(1).R1(1)).^2 + ...
    (trial_points(:,2)./values(1).R2(1)).^2;
vals(vals>1) = NaN;
[~,idx] = max(vals);
dist2 = trial_distances(idx);
ratio = dist1/dist2;

newTraj = zeros(data_size, 3);
newTraj(1,:) = initPoint;

RR_vec = RR_vec/norm(RR_vec);

%Set time-step and flag for debugging
jj2 = 1;
plotting_flag = false; %Put breakpoint at Line 183 before changing this
for jj = jj2+1:jj2:data_size
    
    %Find max distance
    trial_distances = linspace(0, values(1).R1(jj), 1000)';
    trial_points = repmat(RR_vec, length(trial_distances), 1) ...
        .*repmat(trial_distances, 1, 2);
    vals = (trial_points(:,1)./values(1).R1(jj)).^2 + ...
        (trial_points(:,2)./values(1).R2(jj)).^2;
    vals(vals>1) = NaN;
    [~,idx] = max(vals);
    dist2 = trial_distances(idx);
    
    %Get reference frame from previous point and current point
    v1 = values(1).T(jj-jj2,:); w1 = values(1).T(jj,:);
    v2 = values(1).N(jj-jj2,:); w2 = values(1).N(jj,:);
    v3 = values(1).B(jj-jj2,:); w3 = values(1).B(jj,:);
    
    %Get unit vector from the previous mean and previous point in
    %trajectory (for direction)
    vector1 = (newTraj(jj-jj2,:) - mean_vec(jj-jj2,:));
    vector1 = vector1/norm(vector1);
    
    %Move vector to new TNB frame
    vector2 = dot(vector1, v2)*w2 + dot(vector1, v3)*w3;
    
    %Rotates unit vector from previous reference frame to have the same
    %orientation based on new reference frame
  
    %Get the difference in the alphas
    new_alpha = (values(1).alpha(jj) - values(1).alpha(jj-jj2));
    
    %To fix the occasional rotation by 180 degrees
    if new_alpha < -pi/2
        new_alpha = new_alpha + pi;
    elseif new_alpha > pi/2
        new_alpha = new_alpha - pi;
    end
    
    %Rotate around T2 vector by the difference in the alphas
    vector3 = (AxelRot(vector2', new_alpha*180/pi, w1, [0 0 0]))';
    
    %Scales the magnitude to have the same ratio when comparing to the
    %maximum distance
    vector4 = vector3/norm(vector3)*ratio*dist2;
    
    newTraj(jj,:) = vector4 + mean_vec(jj,:);
    if plotting_flag %Used for debugging purposes
        hold on; 
        plot3(mean_vec(jj-jj2,1), mean_vec(jj-jj2,2), mean_vec(jj-jj2,3), 'r+', 'Linewidth', 2);
        x = values(1).R1(jj-jj2)*cos(2*pi*tc)*cos(values(1).alpha(jj-jj2)) - values(1).R2(jj-jj2)*sin(2*pi*tc)*sin(values(1).alpha(jj-jj2));
        y = values(1).R1(jj-jj2)*cos(2*pi*tc)*sin(values(1).alpha(jj-jj2)) + values(1).R2(jj-jj2)*sin(2*pi*tc)*cos(values(1).alpha(jj-jj2));
        C = repmat(mean_vec(jj-jj2,:)',1,length(tc)) + values(1).N(jj-jj2,:)'*x + values(1).B(jj-jj2,:)'*y;
        plot3(C(1,:), C(2,:), C(3,:), 'r');
        plot3(newTraj(jj-jj2, 1), newTraj(jj-jj2,2), newTraj(jj-jj2,3), 'r+', 'Linewidth', 2);
        plot3([mean_vec(jj-jj2,1), newTraj(jj-jj2,1)], [mean_vec(jj-jj2,2), newTraj(jj-jj2,2)], [mean_vec(jj-jj2,3), newTraj(jj-jj2,3)], 'r');

        r1_vec2 = (AxelRot(values(1).N(jj-jj2,:)', values(1).alpha(jj-jj2)*180/pi, values(1).T(jj-jj2,:), [0 0 0]))';
        r1_point2 = r1_vec2/norm(r1_vec2)*values(1).R1(jj-jj2) + mean_vec(jj-jj2,:);
        plot3([mean_vec(jj-jj2,1), r1_point2(1)], [mean_vec(jj-jj2,2), r1_point2(2)], [mean_vec(jj-jj2,3), r1_point2(3)], 'k', 'Linewidth', 1.2);
        r2_vec2 = (AxelRot(values(1).N(jj-jj2,:)', values(1).alpha(jj-jj2)*180/pi+90, values(1).T(jj-jj2,:), [0 0 0]))';
        r2_point2 = r2_vec2/norm(r2_vec2)*values(1).R2(jj-jj2) + mean_vec(jj-jj2,:);
        plot3([mean_vec(jj-jj2,1), r2_point2(1)], [mean_vec(jj-jj2,2), r2_point2(2)], [mean_vec(jj-jj2,3), r2_point2(3)], 'k', 'Linewidth', 1.2);
        r1_vec2 = (AxelRot(values(1).N(jj-jj2,:)', values(1).alpha(jj-jj2)*180/pi + 180, values(1).T(jj-jj2,:), [0 0 0]))';
        r1_point2 = r1_vec2/norm(r1_vec2)*values(1).R1(jj-jj2) + mean_vec(jj-jj2,:);
        plot3([mean_vec(jj-jj2,1), r1_point2(1)], [mean_vec(jj-jj2,2), r1_point2(2)], [mean_vec(jj-jj2,3), r1_point2(3)], 'k');
        r2_vec2 = (AxelRot(values(1).N(jj-jj2,:)', values(1).alpha(jj-jj2)*180/pi+90+180, values(1).T(jj-jj2,:), [0 0 0]))';
        r2_point2 = r2_vec2/norm(r2_vec2)*values(1).R2(jj-jj2) + mean_vec(jj-jj2,:);
        plot3([mean_vec(jj-jj2,1), r2_point2(1)], [mean_vec(jj-jj2,2), r2_point2(2)], [mean_vec(jj-jj2,3), r2_point2(3)], 'k');

        plot3(mean_vec(jj,1), mean_vec(jj,2), mean_vec(jj,3), 'b+', 'Linewidth', 2);
        x = values(1).R1(jj)*cos(2*pi*tc)*cos(values(1).alpha(jj)) - values(1).R2(jj)*sin(2*pi*tc)*sin(values(1).alpha(jj));
        y = values(1).R1(jj)*cos(2*pi*tc)*sin(values(1).alpha(jj)) + values(1).R2(jj)*sin(2*pi*tc)*cos(values(1).alpha(jj));
        C = repmat(mean_vec(jj,:)',1,length(tc)) + values(1).N(jj,:)'*x + values(1).B(jj,:)'*y;
        plot3(C(1,:), C(2,:), C(3,:), 'b');
        plot3(newTraj(jj, 1), newTraj(jj,2), newTraj(jj,3), 'b+', 'Linewidth', 2);
        plot3([mean_vec(jj,1), newTraj(jj,1)], [mean_vec(jj,2), newTraj(jj,2)], [mean_vec(jj,3), newTraj(jj,3)], 'b');

        r1_vec2 = (AxelRot(values(1).N(jj,:)', values(1).alpha(jj)*180/pi, values(1).T(jj,:), [0 0 0]))';
        r1_point2 = r1_vec2/norm(r1_vec2)*values(1).R1(jj) + mean_vec(jj,:);
        plot3([mean_vec(jj,1), r1_point2(1)], [mean_vec(jj,2), r1_point2(2)], [mean_vec(jj,3), r1_point2(3)], 'k', 'Linewidth', 1.2);
        r2_vec2 = (AxelRot(values(1).N(jj,:)', values(1).alpha(jj)*180/pi+90, values(1).T(jj,:), [0 0 0]))';
        r2_point2 = r2_vec2/norm(r2_vec2)*values(1).R2(jj) + mean_vec(jj,:);
        plot3([mean_vec(jj,1), r2_point2(1)], [mean_vec(jj,2), r2_point2(2)], [mean_vec(jj,3), r2_point2(3)], 'k', 'Linewidth', 1.2);
        r1_vec2 = (AxelRot(values(1).N(jj,:)', values(1).alpha(jj)*180/pi + 180, values(1).T(jj,:), [0 0 0]))';
        r1_point2 = r1_vec2/norm(r1_vec2)*values(1).R1(jj) + mean_vec(jj,:);
        plot3([mean_vec(jj,1), r1_point2(1)], [mean_vec(jj,2), r1_point2(2)], [mean_vec(jj,3), r1_point2(3)], 'k');
        r2_vec2 = (AxelRot(values(1).N(jj,:)', values(1).alpha(jj)*180/pi+90+180, values(1).T(jj,:), [0 0 0]))';
        r2_point2 = r2_vec2/norm(r2_vec2)*values(1).R2(jj) + mean_vec(jj,:);
        plot3([mean_vec(jj,1), r2_point2(1)], [mean_vec(jj,2), r2_point2(2)], [mean_vec(jj,3), r2_point2(3)], 'k');
        axis equal;
        hold off;
        clf;
    end
end

end