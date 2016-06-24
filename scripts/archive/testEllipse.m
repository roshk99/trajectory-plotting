clc;clear;close all;
max_point = [-0.262726823499201,0.675400613641133,0.406002945570792;-0.357112776023504,0.0367393149922423,0.0784978028808475;-0.341257126372320,0.837622670019911,0.188414043984817;0.00694650217388605,-0.648276799904241,0.231649620442214];
T = [-0.973564193596502;-0.210611588960103;0.0884054271099588];
point1 = [-0.274496641834378,0.705204680822199,0.352382456396313];
N = [-0.207732348240021;0.977338637531833;0.0406996201715280];
B = [-0.0949738513564651;0.0212590259212638;-0.995252742460627];
d = -dot(point1,T);
[xx,yy] = ndgrid(-2:0.1:2, -2:0.1:2);
z = (-T(1)*xx - T(2)*yy - d)/T(3);

subplot(1,2,1);
plot3(max_point(:,1), max_point(:,2), max_point(:,3), 'r+', 'Linewidth', 4);
hold on;
plot3(point1(1), point1(2), point1(3), 'b+', 'Linewidth', 5);
% mesh(xx,yy,z);
axis equal;
view([0.9 90]);
set(gcf, 'Position', get(0, 'Screensize'));

% uv = null(T.');
% uv = [N B];
% new_points = max_point*uv;
% new_center = point1*uv;

for k=1:size(max_point,1)
new_points(k,:) = [dot(N, max_point(k,:)-point1), dot(B, max_point(k,:)-point1)];
end
new_center = [0,0];

[a1,b1,vect_a,vect_b,alpha] = estimateEllipse(new_points, new_center);
tc = 0:0.01:1; % arc-length on the circle
x = a1*cos(2*pi*tc)*cos(alpha) - b1*sin(2*pi*tc)*sin(alpha);
y = a1*cos(2*pi*tc)*sin(alpha) + b1*sin(2*pi*tc)*cos(alpha);
C = repmat(point1',1,length(tc)) + N*x + B*y;

plot3(C(1,:), C(2,:), C(3,:), 'g', 'Linewidth', 2);
plot3([point1(1) point1(1)+1.5*N(1)], [point1(2) point1(2)+1.5*N(2)], [point1(3) point1(3)+1.5*N(3)], 'c', 'Linewidth', 2);
plot3([point1(1) point1(1)+1.5*B(1)], [point1(2) point1(2)+1.5*B(2)], [point1(3) point1(3)+1.5*B(3)], 'c', 'Linewidth', 2);
hold off;

subplot(1,2,2);
plot(new_points(:,1), new_points(:,2), 'r+', 'Linewidth', 5);
hold on;
plot(new_center(1), new_center(2), 'b+', 'Linewidth', 5);
plot(new_center(1)+x, new_center(2)+y, 'g', 'Linewidth', 5);
hold off;
axis equal;