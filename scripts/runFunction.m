function runFunction(set_num, fit_type, plot_method, how_many, handles)
%Parameters

%If data needs to be shifted to a specific end point
end_point = [Inf,Inf,Inf;0,0,0;Inf,Inf,Inf;0,0,0;Inf,Inf,Inf; Inf Inf Inf];

%Plotting Flags
plotRawTrajectories = false;
plotSmoothTrajectories = false;
plotComponents = false;
plotSurface = true;
plotNewTrajectories = true;
plotBoundaries = false;

%Parsing Parameters
span = 0.15; %used for smoothing
tol = 5e-4; %used for cutting off unnecessary data in beginning and end

%Resampling parameters
point_num = 250;
neighborhood = point_num*0.1; %Number of points to consider on either 
                              %side of a point when finding boundaries
%Number of points to cut off on either side of trajectories
if set_num == 0
    idx1 = 1; idx2 = 0;
else
    idx1 = round(point_num*0.1);
    idx2 = round(point_num*0.1);
end

%Reproduction Parameters
random_flag = true;
trajectory_num = 5;

if set_num > 0
    [raw_data, data] = parseTrajectory(set_num, end_point(set_num+1, :), ...
        tol, span, point_num);
    
else
    load artificial_data;
    plotRawTrajectories = false;
    raw_data = data;
end

[values, canal] = encodeTrajectory(data, how_many, idx1, ...
    idx2, neighborhood, fit_type);

if strcmp(fit_type, 'circles')
    traj = createTrajectory_circle(values, random_flag, trajectory_num);
else
    traj = createTrajectory_ellipse(values, random_flag, trajectory_num);
    values(1).Router = [];
end

mean_vals = {};
mean_vals{1} = values(1).xmean;
mean_vals{2} = values(1).ymean;
mean_vals{3} = values(1).zmean;

%Cut original data to match the rest of the data used for plotting
for ii=1:length(data)
    data{ii} = data{ii}(idx1:end-idx2,:);
end

plotting(raw_data, data, mean_vals, canal, traj, values(1).Router, ...
    plotRawTrajectories, plotSmoothTrajectories, plotComponents, ...
    plotSurface, plotNewTrajectories, plotBoundaries, how_many, ...
    plot_method, handles)
styles(set_num, handles); %Styles the plot (set the title in this function)
end