function runFunction(set_num, fit_type, plot_method, how_many, handles)
    %Parameters
    end_point = [Inf,Inf,Inf;0,0,0;Inf,Inf,Inf;0,0,0;Inf,Inf,Inf; Inf Inf Inf];
    idx1 = [1,200,250,275,225,1];
    idx2 = [0,200,200,200,225,0];
    plotRawTrajectories = false;
    plotSmoothTrajectories = false;
    plotComponents = false;
    plotSurface = true;
    plotNewTrajectories = false;
    point_num = 2000;
    neighborhood = point_num*0.1;
    random_flag = 1;
    span = 0.15;
    tol = 5e-4;
    trajectory_num = 3;
    view_vec = [159.2,26;19.6,-6.8;-122.8,24.4;-48.8,3.6; 14 4.4; 49.2 22];
    
    if set_num > 0
        [raw_data, data] = parseTrajectory(set_num, end_point(set_num+1, :), ...
            tol, span, point_num);
        
    else
        load artificial_data;
        plotRawTrajectories = false;
        raw_data = {};
    end

    [values, canal] = encodeTrajectory(data, how_many, idx1(set_num+1), ...
        idx2(set_num+1), neighborhood, fit_type);
    
    if strcmp(fit_type, 'circles')
        traj = createTrajectory_circle(values, random_flag, trajectory_num);
    else
        traj = {};
    end
    
    mean_vals = {};
    mean_vals{1} = values(1).xmean;
    mean_vals{2} = values(1).ymean;
    mean_vals{3} = values(1).zmean;
    
    for ii=1:length(data)
        data{ii} = data{ii}(idx1(set_num+1):end-idx2(set_num+1),:);
    end
    plotting(raw_data, data, mean_vals, canal, traj, values(1).Router, ...
        plotRawTrajectories, plotSmoothTrajectories, plotComponents, ...
        plotSurface, plotNewTrajectories, how_many, plot_method, ...
        set_num, fit_type, view_vec(set_num+1,:), handles);
end