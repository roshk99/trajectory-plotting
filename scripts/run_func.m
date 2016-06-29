function traj = run_func(set_num, fit_type, plot_method, how_many, handles)
    %Parameters
    end_point = [Inf,Inf,Inf;0,0,0;Inf,Inf,Inf;0,0,0;Inf,Inf,Inf];
    idx1 = [1,200,250,275,200];
    idx2 = [0,200,200,200,200];
    plotNewTrajectories = false;
    plotRawTrajectories = true;
    plotSurface = false;
    point_num = 5000;
    random_flag = 1;
    span = 0.15;
    tol = 5e-4;
    trajectories = 1;
    view_vec = [-42,-18.8;19.6,-6.8;-91.2,-3.6;-48.8,3.6; -57.6 23.6];
    
    if set_num > 0
        data = parsetrajectory(set_num, end_point(set_num+1, :), ...
            tol, span, point_num, plotRawTrajectories);
    else
        load artificial_data;
    end
    
    traj = create_trajectory(data, idx1(set_num+1), ...
        idx2(set_num+1), set_num, how_many, ...
        view_vec(set_num+1,:), random_flag, ...
        trajectories, plotSurface, plotNewTrajectories, fit_type, ...
        plot_method, handles);
end