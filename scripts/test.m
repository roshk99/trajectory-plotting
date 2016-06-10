mdl_jaco;
path = [0.0690570076448947,0.117485534182401,0.295923229572532];
offset = [0.25 0.25 0.25];
path = path + offset;
figure;
plot3(path(1), path(2), path(3), 'k+', 'LineWidth', 8);
hold on;
q0 = [123 -173 81.3 0 0 0]*pi/180;
% jaco.teach();
axis([-0.5 2 -1 1 -0.5 1])
Tp = transl(path);
%q = jaco.ikine(Tp, 'tol', 1e-2, 'alpha', 0.01, 'ilimit', 3000)
q = jaco.ikcon(Tp, q0);
jaco.plot(q);
hold off;
%   Columns 1 through 5
% 
%   182.9913 -340.7323  -65.2656   31.6502 -101.9120
% 
%   Column 6
% 
%   108.6614