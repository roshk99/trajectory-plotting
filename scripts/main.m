function main()
    %Global parameters
    point_num = 1000; %number points desired in final trajectory
    span = 0.15; %affects the smoothness of the data
    tol = 0.0005; %the lowest variation desired in initial and final points
    plot_bool = false;
    record_bool = true;
    
    %Dataset1
    set_num1 = 1;
    end_point1 = [0 0 0]; %point to shift the end of all trajectories
    offset_vec1 = 0.25*ones(4, 3); %shift the trajectories
    object_placement1 = [0.25 0.25 0.25];
    axis1 = [-0.5 0.5 -0.5 0.5 0 0.8];
    set1_title = 'Robot approaching the ball with varying trajectories';
    smooth_data1 = parsetrajectory(set_num1, end_point1, tol, span, ...
        point_num, plot_bool);
    animate(set_num1, smooth_data1, offset_vec1, axis1, ...
        object_placement1, record_bool, set1_title);
    
    %Dataset2
    set_num2 = 2;
    end_point2 = false;
    offset_vec2 = [0 0.05 -0.25;...
                  0 0.05 -0.25;...
                  0.1 0 -0.29;...
                  0 -0.05 -0.25];
    object_placement2 = [-0.4 0.11 -0.15];
    axis2 = [-1 1 -1 1 -0.5 0.5];
    set2_title = 'Robot picking up the ball with varying trajectories';
    smooth_data2 = parsetrajectory(set_num2, end_point2, tol, span, ...
        point_num, plot_bool);
    animate(set_num2, smooth_data2, offset_vec2, axis2, ...
        object_placement2, record_bool, set2_title);
end