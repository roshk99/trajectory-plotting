function main()
    %Global parameters
    point_num = 1000; %number points desired in final trajectory
    span = 0.15; %affects the smoothness of the data
    tol = 0.0005; %the lowest variation desired in initial and final points
    plot_bool = false;
    record_bool = false;
    q0 = [123 -173 81.3 90 -90 0]*pi/180;
    
    %Dataset1
    set_num1 = 1;
    end_point1 = [0 0 0]; %point to shift the end of all trajectories
    offset_vec1 = [0.5 0.25 0.25; ...
                    0.5 0.25 0.25;...
                    0.5 0.25 0.25;...
                    0.5 0.25 0.25];
    object_placement1 = [0.5 0.25 0.25];
    view1 = [-127.5, 30];
    canalview1 = [19.6 -6.8];
    axis1 = [-0.25 1 -0.5 0.6 -0.5 0.8];
    set1_title = 'Robot approaching the ball with varying trajectories';
    smooth_data1 = parsetrajectory(set_num1, end_point1, tol, span, ...
        point_num, plot_bool);
    plot_boundaries(smooth_data1, set_num1, canalview1);
    animate(set_num1, smooth_data1, offset_vec1, view1, axis1, ...
        object_placement1, record_bool, set1_title, q0);
    
    %Dataset2
    set_num2 = 2;
    end_point2 = false;
    offset_vec2 = [-0.1 0.05 -0.25;...
              -0.1 0.05 -0.25;...
              0 0 -0.29;...
              -0.1 -0.05 -0.25];
    view2 = [-127.5, 30];
    canalview2 = [0 90];
    object_placement2 = [-0.5 0.03 -0.15];
    axis2 = [-1 1 -1 1 -0.5 0.5];
    set2_title = 'Robot picking up the ball with varying trajectories';
    smooth_data2 = parsetrajectory(set_num2, end_point2, tol, span, ...
        point_num, plot_bool);
    plot_boundaries(smooth_data2, set_num2, canalview2);
    animate(set_num2, smooth_data2, offset_vec2, view2, axis2, ...
        object_placement2, record_bool, set2_title, q0);

    %Dataset3
    set_num3 = 3;
    end_point3 = [0 0 0];
    offset_vec3 = [0.25 0.25 0.25;
                   0.25 0.25 0.25;
                   0.25 0.25 0.25;
                   0.25 0.25 0.25;
                   0.25 0.25 0.25];
    object_placement3 = [0.6 -0.3 0.7];
    view3 = [-200 30];
    canalview3 = [0 90];
    axis3 = [-0.5 0.75 -0.9 0.5 -1 0.9];
    q0 = [270 180 90 0 90 45]*pi/180;
    set3_title = 'Robot moving between two points';
    smooth_data3 = parsetrajectory(set_num3, end_point3, tol, span, ...
        point_num, plot_bool);
    plot_boundaries(smooth_data3, set_num3, canalview3);
    animate(set_num3, smooth_data3, offset_vec3, view3, axis3, ...
        object_placement3, record_bool, set3_title, q0);
end