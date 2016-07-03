function runFunction(set_num, fit_type, plot_method, how_many, handles)
    %Parameters
    end_point = [Inf,Inf,Inf;0,0,0;Inf,Inf,Inf;0,0,0;Inf,Inf,Inf];
    idx1 = [1,200,250,275,225];
    idx2 = [0,200,200,200,225];
    plotRawTrajectories = false;
    plotSmoothTrajectories = false;
    plotComponents = false;
    plotSurface = true;
    plotNewTrajectories = false;
    point_num = 5000;
    random_flag = 1;
    span = 0.15;
    tol = 5e-4;
    trajectory_num = 1;
    view_vec = [-42,-18.8;19.6,-6.8;-91.2,-3.6;-48.8,3.6; -20.4 0.4];
    
    if set_num > 0
        [raw_data, data] = parseTrajectory(set_num, end_point(set_num+1, :), ...
            tol, span, point_num);
        
    else
        load artificial_data;
        plotRawTrajectories = false;
        raw_data = {};
    end
    
    [values, canal] = encodeTrajectory(data, how_many, idx1(set_num+1), ...
        idx2(set_num+1), fit_type);
    
    if strcmp(fit_type, 'circles')
        traj = createTrajectory_circle(values, random_flag, trajectory_num);
    end
    
    mean_vals = {};
    mean_vals{1} = values(1).xmean;
    mean_vals{2} = values(1).ymean;
    mean_vals{3} = values(1).zmean;
    
    plotting(raw_data, data, mean_vals, canal, traj, plotRawTrajectories, ...
        plotSmoothTrajectories, plotComponents, plotSurface, ...
        plotNewTrajectories, how_many, plot_method, set_num, fit_type, ...
        view_vec(set_num+1,:), handles);   
end