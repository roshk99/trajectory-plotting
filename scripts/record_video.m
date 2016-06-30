function record_video(data, view_vec)

offset_vec = [0 0 0; ...
              0 0 0;...
              0 0 0;...
              0 0 0;...
              0 0 0];


%Plot object
C = [-0.23,-0.35,0.53];
[x,y,z] = sphere(20);
scale = 0.18;
surf(scale*x+C(1,1),scale*y+C(1,2),scale*z+C(1,3));
hold on;
view(view_vec);
grid on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
set(gcf, 'Position', get(0, 'Screensize'));

cc = jet(length(data));

%Animate first two trajectories
jaco = mdl_jaco_custom();
for ind1 = 1:2
    path = data{ind1};
    
    path = [path(:,1)+offset_vec(ind1, 1)...
        path(:,2)+offset_vec(ind1, 2)...
        path(:,3)+offset_vec(ind1, 3)];
    
    plot3(path(:,1), path(:,2), path(:,3), 'color', ...
        cc(ind1, :), 'LineWidth', 2);
    
    path = path(1:30:end, :);
    Dt = 1;
    Qdmax = [0.5 0.5 0.3];
    Tacc = 0.01;
    p0xyz = [0 0 0];
    p = mstraj(path, Qdmax, [], p0xyz, Dt, Tacc);
    
    Tp = transl(p);
    q = jaco.ikunc(Tp);
    directoryname = sprintf('set%i - run%i', set_num, ind1);
    jaco.plot(q);
end
hold off;

end