function canal_visualization(set_to_run, fit_type, how_many, ...
    boundary_method, cross_section_method, plot_method, handles)
% -----------------------------------------------------------------------
% A function that generates a canal surface for a dataset
%
% Inputs:
%   set_to_run: number of dataset to use
%   fit_type: 'circle', 'ellipse', or 'bspline' based on type of surface
%             desired
%   boundary_method: 1 or 2 depending on which method is used to calculate
%                    the boundary (only for circles)
%   cross_section_method: 1 or 2 depending on whether TNB frames are used
%                         to calculate cross-section alignment or the
%                         cross-section itself
%   plot_method: 'circles' or 'surface' based on desired plot
%   handles: the gui object for canal visualization gui
%
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------

load global_parameters

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
        plotNewTrajectories, fit_type, boundary_method, ...
        cross_section_method, plot_method, handles);    
elseif set_to_run == 1
    %Dataset1    
    smooth_data1 = parsetrajectory(1, [0 0 0], tol, span, ...
        point_num, plotRawTrajectories);
    traj1 = create_trajectory(smooth_data1, 100, 100, 1, ...
        how_many, [19.6 -6.8], true, 10, plotSurface, ...
        plotNewTrajectories, fit_type, boundary_method, ...
        cross_section_method, plot_method, handles);  
elseif set_to_run == 2
    %Dataset2
    smooth_data2 = parsetrajectory(2, false, tol, span, ...
        point_num, plotRawTrajectories);
    traj2 = create_trajectory(smooth_data2, 100, 100, 2, ...
        how_many, [0 90], true, 10, plotSurface, ...
        plotNewTrajectories, fit_type, boundary_method, ...
        cross_section_method, plot_method, handles);  
elseif set_to_run == 3
    %Dataset3
    smooth_data3 = parsetrajectory(3, [0 0 0], tol, span, ...
        point_num, plotRawTrajectories);
    traj3 = create_trajectory(smooth_data3, 100, 100, 3, ...
        how_many, [-48.8 3.6], true, 10, plotSurface, ...
        plotNewTrajectories, fit_type, boundary_method, ...
        cross_section_method, plot_method, handles);  
end
end