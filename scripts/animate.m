

function animate(set_num, data, offset_vec, axis_vec, ...
    object_placement, record_bool, plot_title)
    mdl_puma560;
    
    %Plot the object
    r = 0.06;
    [x,y,z] = sphere(50);
    x0 = object_placement(1); y0 = object_placement(2); 
    z0 = object_placement(3);
    x = x*r + x0; y = y*r + y0; z = z*r + z0 - 1.5*r;
    surface(x,y,z);
    hold on;
    
    cc = jet(length(data));
    for ind1 = 1:length(data)
        path = data{ind1};

        path = [path(:,1)+offset_vec(ind1, 1)...
                path(:,2)+offset_vec(ind1, 2)...
                path(:,3)+offset_vec(ind1, 3)];

        plot3(path(:,1), path(:,2), path(:,3), 'color', ...
            cc(ind1, :), 'LineWidth', 2);
        hold on;
        axis(axis_vec);
        view([-127.5, 30]);
        grid on;
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        title(plot_title);
        
        path = path(1:10:end, :);
        Dt = 1;
        Qdmax = [0.5 0.5 0.3];
        Tacc = 0.01;
        p0xyz = [0 0 0];
        p = mstraj(path, Qdmax, [], p0xyz, Dt, Tacc);

        Tp = transl(p);
        q = p560.ikine6s(Tp, 'run');
        directoryname = sprintf('set%i - run%i', set_num, ind1);
        if record_bool
            p560.plot(q, 'movie', directoryname);
        else
            p560.plot(q);
        end
    end
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title(plot_title);
end