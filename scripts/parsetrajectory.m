function smooth_data = parsetrajectory(set_num, end_point, tol, ...
    span, point_num, plot_bool)   
    %Get the x,y,z trajectories
    data = get_raw_trajectories(set_num);
    
    %Cut data
    cut_data = get_cut_trajectories(data, end_point, tol);

    %Fit Data
    fit_data = get_fit_trajectories(cut_data, point_num);

    %Smooth data
    smooth_data = get_smooth_trajectories(fit_data, span);
                  
    %Plots
    if plot_bool
        plot_data(fit_data, sprintf('Raw Data for Set %i', set_num));
        plot_data(smooth_data, sprintf('Smoothed Data for Set %i', set_num));
    end
end

function data_vars = get_raw_trajectories(set_num)
    %Load variables
    if set_num == 1
        load dataset1;
        raw_vars = {dataset1_1, dataset1_3, dataset1_6, dataset1_7};
    elseif set_num == 2
        load dataset2;
        raw_vars = {dataset2_4, dataset2_5, dataset2_6, dataset2_7};
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
        for ind2 = 1:3
            traj = traj(abs(traj(:,ind2) - ...
                traj(1,ind2)*ones(length(traj(:,ind2)), 1)) > ...
                tol*ones(length(traj(:,ind2)), 1), :);
            traj = traj(abs(traj(:,ind2) - ...
                traj(end,ind2)*ones(length(traj(:,ind2)), 1)) > ...
                tol*ones(length(traj(:,ind2)), 1), :);
        end
        
        %Shift data
        if length(end_point) == 3
            cut_data{ind1} = [traj(:,1) - traj(end, 1) + end_point(1),...
                traj(:,2) - traj(end, 2) + end_point(2),...
                traj(:,3) - traj(end, 3) + end_point(3)];
        else
            cut_data{ind1} = traj;
        end
    end
end

function fit_data = get_fit_trajectories(data, point_num)
    
    new_t_vec = (0:point_num-1)';
    
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

function plot_data(data, plot_title)
    figure;
    cc = jet(length(data));
    if isempty(plot_title)
        subplot(3,1,1); hold all; title('X Component');
    else
        subplot(3,1,1); hold all; 
        title(sprintf('%s - X Component', plot_title));
    end
    subplot(3,1,2); hold all; title('Y Component');
    subplot(3,1,3); hold all; title('Z Component');
    
    legend_vec = {};
    for ind1 = 1:length(data)
        tvec = -(length(data{ind1}(:,1))-1):0;
        subplot(3,1,1);
        plot(tvec, data{ind1}(:,1), 'color', cc(ind1,:), ...
            'linewidth', 1.2);
        subplot(3,1,2);
        plot(tvec, data{ind1}(:,2), 'color', cc(ind1,:), ...
            'linewidth', 1.2);
        subplot(3,1,3);
        plot(tvec, data{ind1}(:,3), 'color', cc(ind1,:), ...
            'linewidth', 1.2);
        legend_vec{ind1} = sprintf('Series %i', ind1);
    end
    legend(legend_vec);
    subplot(3,1,1); hold off; axis tight;
    subplot(3,1,2); hold off; axis tight;
    subplot(3,1,3); hold off; axis tight;
    
    figure;
    for ind1 = 1:length(data)
       plot3(data{ind1}(:, 1), data{ind1}(:,2), data{ind1}(:,3), ...
           'color', cc(ind1,:), 'linewidth', 1.5);
       hold on;
    end
    hold off;
    if isempty(plot_title)
        title('3D Trajectories');
    else
        title(sprintf('%s - 3D Trajectories', plot_title));
    end
    legend(legend_vec);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on;
    axis equal;
    
end