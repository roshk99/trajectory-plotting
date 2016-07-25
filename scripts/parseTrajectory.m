function [cut_data, smooth_data] = parseTrajectory(set_num, end_point, tol, ...
    span, point_num)
% -----------------------------------------------------------------------
% A function to analyze raw data and output smooth data
%
% Inputs:
%   set_num: integer used to label the plot
%   end_point: point to which all trajectories will be shifted, false if
%              shifting not desired
%   tol: amount of variation desired in beginning and end of trajectories
%   span: input to the smoothing function, degree of smoothness
%   point_num: number of points desired in each smoothed trajectories,
%              trajectories will be resampled to achieve this
%
% Output:
%   cut_data: same as smooth_data without fitting and smoothing
%   smooth_data: cell with each element point_numx3 vector of smoothed data
%
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------

%Get the x,y,z trajectories
data = get_raw_trajectories(set_num);

%Cut data
cut_data = get_cut_trajectories(data, end_point, tol);

%Fit Data
fit_data = get_fit_trajectories(cut_data, point_num);

%Smooth data
smooth_data = get_smooth_trajectories(fit_data, span);
end

function data_vars = get_raw_trajectories(set_num)
%Load variables
if set_num == 1
    load dataset1;
    raw_vars = {dataset1_1, dataset1_3, dataset1_6, dataset1_7};
elseif set_num == 2
    load dataset2;
    raw_vars = {dataset2_4, dataset2_5, dataset2_6, dataset2_7};
elseif set_num == 3
    load dataset3;
    raw_vars = {dataset3_1, dataset3_4, dataset3_5, dataset3_6, dataset3_7};
elseif set_num == 4
    load dataset4;
    raw_vars = {dataset4_3, dataset4_4, dataset4_5, dataset4_6, dataset4_7};
elseif set_num == 5
    load dataset5;
    raw_vars = {dataset5_1, dataset5_2, dataset5_3, dataset5_4};
else
    fprintf('Error Occurred - Dataset %i not available', set_num);
    return;
end

%Initialize vectors
data_vars = {};

%Store the x,y,z position for each trajectory
for ind1 = 1:length(raw_vars)
    data_vars{ind1} = raw_vars{ind1}(:, 8:10);
end
end

function cut_data = get_cut_trajectories(data, end_point, tol)
cut_data = {};
arc_length_vec = zeros(length(data), 1);
for ind1 = 1:length(data)
    traj = data{ind1};
    
    %Eliminate nonvarying poritions in beginning and end
%     for ind2 = 1:3
%         traj = traj(abs(traj(:,ind2) - ...
%             traj(1,ind2)*ones(length(traj(:,ind2)), 1)) > ...
%             tol*ones(length(traj(:,ind2)), 1), :);
%         traj = traj(abs(traj(:,ind2) - ...
%             traj(end,ind2)*ones(length(traj(:,ind2)), 1)) > ...
%             tol*ones(length(traj(:,ind2)), 1), :);
%     end

    traj = traj(sum(abs(traj - repmat(traj(1,:), length(traj), 1)) ...
        > repmat(tol, length(traj), 3), 2) > 0, :);
    traj = traj(sum(abs(traj - repmat(traj(end,:), length(traj), 1)) ...
        > repmat(tol, length(traj), 3), 2) > 0, :);


    %Shift data
    if end_point ~= [Inf, Inf, Inf]
        cut_data{ind1} = [traj(:,1) - traj(end, 1) + end_point(1),...
            traj(:,2) - traj(end, 2) + end_point(2),...
            traj(:,3) - traj(end, 3) + end_point(3)];
    else
        cut_data{ind1} = traj;
    end
end
end

function fit_data = get_fit_trajectories(data, point_num)

%Fit data and get new warped data
for ind1 = 1:length(data)
    traj_vec = data{ind1};
    t_vec = (0:length(traj_vec)-1)';
    new_t_vec = linspace(0, length(traj_vec)-1, point_num)';
    
    fit_data_vec = zeros(point_num, 3);
    for ind2 = 1:3
        fitobject = fit(t_vec, traj_vec(:,ind2), 'smoothingspline');
        new_traj = fitobject(new_t_vec);
        fit_data_vec(:, ind2) = new_traj;
    end
    fit_data{ind1} = fit_data_vec;
end

end

function smooth_data = get_smooth_trajectories(data, span)
smooth_data = {};

for ind1 = 1:length(data)
    traj = data{ind1};
    
    smooth_traj = zeros(size(traj));
    for ind2 = 1:3
        smooth_traj(:,ind2) = smooth(traj(:,ind2), span);
    end
    smooth_data{ind1} = smooth_traj;
end
end

