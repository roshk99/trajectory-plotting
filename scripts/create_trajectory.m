function trajectories = create_trajectory(data, idx1, idx2, set_num, ...
    howmany, view_vec, randomflag, points, plotSurface, ...
    plotTrajectories, fit_type, canal_method)
% -----------------------------------------------------------------------
% A function that takes smoothed data and generates new trajectories based
% on a canal surface approach
%
% Inputs:
%   data: cell with each element point_numx3 vector of smoothed data
%         (output of parsetrajectory.m)
%   idx1, idx2: values in the beginning and end of trajectories will be
%               eliminated when plotting the surface; data used will be
%               from idx1 to end-idx2
%   set_num: number of dataset, used for plotting
%   howmany: which trajectories to include in analysis (a list of indices)
%   view_vec: specify the view in the 3D plot (a vector with two numbers)
%   randomflag: boolean depending on whether you want random start
%               points for trajectories to be generated
%   points: if randomflag is true, points is an integer with the number of
%           points to generate. if randomflag is false, points is a matrix
%           of dimension point_num x 3 that has the initial points for the
%           new trajectories
%   plotSurface: boolean for plotting the canal surface
%   plotTrajectories: boolean for plotting the generated trajectories
%   fit_type: 'circle', 'ellipse', or 'bspline' based on type of surface
%             desired
%   canal_method: 1 or 2 depending on which method is used to create the
%                 canal
%
% Output:
%   trajectories: cell with each element point_numx3 vector of generated
%                 trajectories
%
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------

%Calculate the boundary (see boundary_calculation function for contents of
%values cell)
values = boundary_calculation(data, set_num, view_vec, howmany, ...
    idx1, idx2, plotSurface, fit_type, canal_method);

%Get the initial mean point
mean_vec = [values(1).xmean, values(1).ymean, values(1).zmean];

%The number of trajectories
data_size = size(values(1).N2, 1);

trajectories = {};
%if randomly generated initial points are desired
if randomflag
    
    %Generates points number of random initial points within the boundary
    %of the canal
    initPoints = zeros(points, 3);
    for ii = 1:points
        if strcmp(fit_type, 'circle')
            initPoints(ii, :) = mean_vec(1,:) + ...
                (2^-.5)*values(1).Router(1)*(-1+2*rand(1))*values(1).N2(1,:) +  ...
                (2^-.5)*values(1).Router(1)*(-1+2*rand(1))*values(1).B2(1,:);
        end
    end
    
    %If initial points are specified
else
    initPoints = points;
end

%Compute the trajectory based on the output from the canal and initial
%points
for ii = 1:length(initPoints)
    if strcmp(fit_type, 'circle')
        trajectories{ii} = compute_trajectory(values, data_size, ...
            initPoints(ii,:));
    end
end

%Plot the new trajectories if desired
if plotTrajectories
    plot_trajectories(trajectories);
end

end

function newTraj = compute_trajectory(values, data_size, initPoint)

%Gets the initial mean point
mean_vec = [values(1).xmean, values(1).ymean, values(1).zmean];

%Gets the distance between the initial point and the mean point
d = norm(mean_vec(1,:) - initPoint);

%Gets the ratio between the distance from the mean to the initial point to
%the outer radius
ratio = d/values(1).Router(1);

newTraj = zeros(data_size, 3);
newTraj(1,:) = initPoint;

%For each point after the initial point
for ii=2:data_size
    
    %Get reference frame from previous point and current point
    v1 = values(1).T(ii-1,:); v1_2 = values(1).T(ii,:);
    v2 = values(1).B2(ii-1,:); v2_2 = values(1).B2(ii,:);
    v3 = values(1).N2(ii-1,:); v3_2 = values(1).N2(ii,:);
    
    %Get unit vector from the previous mean and previous point in
    %trajectory (for direction)
    tmp_unit_vector = (newTraj(ii-1,:) - mean_vec(ii-1,:))...
        /norm((newTraj(ii-1,:) - mean_vec(ii-1,:)));
    
    %Get direction angles between the old and new reference frames
    angles = [dot(v1, v1_2), dot(v2, v2_2), dot(v3, v3_2)];
    
    %Rotates unit vector from previous reference frame to have the same
    %orientation based on new reference frame
    new_vector = tmp_unit_vector.*cos(angles);
    
    %Scales the magnitude to have the same ratio when comparing to new
    %outer radius
    new_vector = new_vector/norm(new_vector)*...
        ratio*values(1).Router(ii);
    
    %Get coordinates of the new point and store it
    newPoint = new_vector + mean_vec(ii,:);
    newTraj(ii,:) = newPoint;
end
end


function plot_trajectories(trajs)
hold on;
for ii = 1:length(trajs)
    plot3(trajs{ii}(:,1), trajs{ii}(:,2), trajs{ii}(:,3), 'g', ...
        'Linewidth', 2);
end
hold off;
set(gcf, 'Position', get(0, 'Screensize'));
end