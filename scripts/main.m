function main()
% -----------------------------------------------------------------------
% A function that analyzes a specific dataset in the specified manner, see
% the parameters below
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------
close all;

set_to_run = 3; %0-3
fit_type = 'circles'; %circles or ellipses
plot_method = 'surface'; %circles or surface
how_many = 1:5;

%Global parameters
load global_parameters;
point_num = 5000;
if set_to_run == 0 %Artificial data
    % Demo-1
    t = linspace(-2,2,100)'; x1 = t; y1 = 1 + x1.^2; z1 = 1*ones(size(t));
    % Demo-2
    x2 = t; y2 = 1 - x2.^2; z2 = 1*ones(size(t));
    % Demo-3
    x3 = t; z3 = 1 + x3.^2; y3 = 1*ones(size(t));
    % Demo-4
    x4 = t; z4 = 1 - x4.^2; y4 = 1*ones(size(t));
    
    data = {[x1 y1 z1], [x2 y2 z2], [x3 y3 z3], [x4 y4 z4]};
    
    traj1 = create_trajectory(data, 1, 0, 0, ...
        how_many, [-42.0000 -18.8000], true, 10, plotSurface, ...
        plotNewTrajectories, fit_type, plot_method, false);
    
elseif set_to_run == 1
%     offset_vec1 = [0.5 0.25 0.25; ...
%         0.5 0.25 0.25;...
%         0.5 0.25 0.25;...
%         0.5 0.25 0.25];
    set1_title = 'Robot approaching the ball with varying trajectories';
    
    smooth_data1 = parsetrajectory(1, [0 0 0], tol, span, ...
        point_num, plotRawTrajectories);
    traj1 = create_trajectory(smooth_data1, 200, 200, 1, ...
        how_many, [19.6 -6.8], true, 10, plotSurface, ...
        plotNewTrajectories, fit_type, plot_method, false);
    %     animate(1, smooth_data1, offset_vec1, [-127.5, 30], ...
    %         [-0.25 1 -0.5 0.6 -0.5 0.8], [0.5 0.25 0.25], record_bool, ...
    %         set1_title, q0);
    %
elseif set_to_run == 2
%     offset_vec2 = [-0.1 0.05 -0.25;...
%         -0.1 0.05 -0.25;...
%         0 0 -0.29;...
%         -0.1 -0.05 -0.25];
    set2_title = 'Robot picking up the ball with varying trajectories';
    smooth_data2 = parsetrajectory(2, false, tol, span, ...
        point_num, plotRawTrajectories);
    traj2 = create_trajectory(smooth_data2, 200, 200, 2, ...
        how_many, [-91.2 -3.6], true, 10, plotSurface, ...
        plotNewTrajectories, fit_type, plot_method, false);
    %     animate(2, smooth_data2, offset_vec2, [-127.5, 30], ...
    %         [-1 1 -1 1 -0.5 0.5], [-0.5 0.03 -0.15], record_bool, ...
    %         set2_title, q0);
elseif set_to_run == 3
%     offset_vec3 = [0.25 0.25 0.25;
%         0.25 0.25 0.25;
%         0.25 0.25 0.25;
%         0.25 0.25 0.25;
%         0.25 0.25 0.25];
    set3_title = 'Robot moving between two points';
    smooth_data3 = parsetrajectory(3, [0 0 0], tol, span, ...
        point_num, plotRawTrajectories);
    traj3 = create_trajectory(smooth_data3, 200, 200, 3, ...
        how_many, [-48.8 3.6], true, 10, plotSurface, ...
        plotNewTrajectories, fit_type, plot_method, false);
    %     animate(3, smooth_data3, offset_vec3, [-200 30], ...
    %         [-0.5 0.75 -0.9 0.5 -1 0.9], [0.6 -0.3 0.7], ...
    %         record_bool, set3_title, q0);
end
end