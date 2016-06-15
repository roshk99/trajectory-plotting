function trajectories = create_trajectory(data, idx1, idx2, set_num, ...
    howmany, view_vec, randomflag, points, plotSurface, plotTrajectories)
    values = boundary_calculation(data, set_num, view_vec, howmany, ...
        idx1, idx2, plotSurface);
    
    mean_vec = [values(1).xmean, values(1).ymean, values(1).zmean];
    data_size = size(values(1).N2, 1);
    trajectories = {};
    if randomflag
        initPoints = zeros(points, 3);
        for ii = 1:points
            initPoints(ii, :) = mean_vec(1,:) + ...
            (2^-.5)*values(1).Router(1)*(-1+2*rand(1))*values(1).N2(1,:) +  ...
            (2^-.5)*values(1).Router(1)*(-1+2*rand(1))*values(1).B2(1,:);
        end
    else
        initPoints = points;
    end
    
    for ii = 1:length(initPoints)
       trajectories{ii} = compute_trajectory(values, data_size, ...
           initPoints(ii,:)); 
    end
    
    if plotTrajectories
       plot_trajectories(trajectories); 
    end
end

function newTraj = compute_trajectory(values, data_size, initPoint) 
    mean_vec = [values(1).xmean, values(1).ymean, values(1).zmean];
    d = norm(mean_vec(1,:) - initPoint);
    ratio = d/values(1).Router(1);

    newTraj = zeros(data_size, 3);
    newTraj(1,:) = initPoint;
    for ii=2:data_size

        v1 = values(1).T(ii-1,:); v1_2 = values(1).T(ii,:);
        v2 = values(1).B2(ii-1,:); v2_2 = values(1).B2(ii,:);
        v3 = values(1).N2(ii-1,:); v3_2 = values(1).N2(ii,:);

        tmp_unit_vector = (newTraj(ii-1,:) - mean_vec(ii-1,:))...
            /norm((newTraj(ii-1,:) - mean_vec(ii-1,:)));
        angles = [dot(v1, v1_2), dot(v2, v2_2), dot(v3, v3_2)];
        new_vector = tmp_unit_vector.*cos(angles);
        new_vector = new_vector/norm(new_vector)*...
            ratio*values(1).Router(ii);
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