function styles(setnum, handles)
if isa(handles, 'struct')
    axes(handles.axes1);
end

if setnum == 1
    %Set 1
    hold on;
    C = [-9.87400965552011e-05;-0.000256255345868566;-0.00211145694817037];
    [x,y,z] = sphere(20);
    scale = 0.008;
    object = surf(scale*x+C(1),scale*y+C(2),scale*z+C(3));
    set(object,'FaceColor','b');
    view([19.6,-6.8]);
elseif setnum == 2
    %Set 2
    view([-296.8,16.4]);
elseif setnum == 4
    %Set 4
    hold on;
    C = [-0.23,-0.35,0.53];
    [x,y,z] = sphere(20);
    scale = 0.17;
    object = surf(scale*x+C(1,1),scale*y+C(1,2),scale*z+C(1,3));
    set(object,'FaceColor','b');
    view([-48.4 9.2]);
end
%All
axis equal; grid on; box on; xlabel('X'); ylabel('Y'); zlabel('Z');
hold off;
title('');
sdf('custom1');

if ~isa(handles, 'struct')
    set(gcf, 'Position', get(0, 'Screensize'));
end
end