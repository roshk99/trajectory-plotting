function values = boundary_calculation(data, set_num, ...
    canalview, howmany, idx1, idx2, plotSurface, fit_type, plot_method, ...
    handles)
% -----------------------------------------------------------------------
% A function that takes smoothed data and generates a canal surface
%
% Inputs:
%   data: cell with each element point_numx3 vector of smoothed data
%         (output of parsetrajectory.m)
%   set_num: number of dataset, used for plotting
%   canalview: specify the view in the 3D plot for the canal (a vector
%              with two numbers
%   howmany: which trajectories to include in analysis (a list of indices)
%   idx1, idx2: values in the beginning and end of trajectories will be
%               eliminated when plotting the surface; data used will be
%               from idx1 to end-idx2
%   plotSurface: boolean for plotting the canal surface
%   fit_type: 'circle', 'ellipse', or 'bspline' based on type of surface
%             desired
%   plot_method: 'circles' or 'surface' based on desired plot
%   handles: the gui object if using the canal visualization gui
%
% Output:
%   values: cell with various canal surface parameters
%           values(1).xmean = mean in x direction (point_num x 1)
%           values(1).ymean = mean in y direction (point_num x 1)
%           values(1).zmean = mean in z direction (point_num x 1)
%           values(1).N = normal component of TNB frame (3 x point_num)
%           values(1).B = binormal component of TNB frame (3 x point_num)
%           values(1).T = tangent component of TNB frame (3 x point_num)
%
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------

%Initialize flags for internal plotting and filtering
plotBoundaries = false;

%Initialize data sizes
numDemos = numel(howmany);
nPoints = size(data{1}, 1);

%Get the X, Y, Z matrices (number of points x number of trajectories)
allXs = zeros(nPoints,numDemos);
allYs = zeros(nPoints,numDemos);
allZs = zeros(nPoints,numDemos);
for ii = 1:numDemos
    allXs(:,ii) = data{ii}(:,1);
    allYs(:,ii) = data{ii}(:,2);
    allZs(:,ii) = data{ii}(:,3);
end
t = linspace(0,1,nPoints).';

%Calculates the mean trajectory
xyz_mean = [mean(allXs, 2), mean(allYs, 2), mean(allZs, 2)];

%Calculate the TNB Frames
[T, N, B] = calculateTNB(t, xyz_mean);

%Finds the boundary values (radii and orientation vectors)
if strcmp(fit_type, 'circles')
    [Router, xyz_distance] = find_boundaries(xyz_mean, allXs, allYs, ...
        allZs);
    
    if plotBoundaries
        %plot_boundaries(data, B, xyz_distance, numDemos, threshold);
    end
elseif strcmp(fit_type, 'ellipses')
    [R1, R2, alpha] = find_boundaries_ellipse(xyz_mean, ...
        allXs, allYs, allZs, T, N);
end


%Get only subset of values within specified boundary
tt = t(idx1:end-idx2);
xx2 = xyz_mean(idx1:end-idx2, 1);
yy2 = xyz_mean(idx1:end-idx2, 2);
zz2 = xyz_mean(idx1:end-idx2, 3);
N = N(:, idx1:end-idx2);
B = B(:, idx1:end-idx2);
values = struct([]);

%Get the canal surface depending on the fit type
if strcmp(fit_type, 'circles')
    RRouter = Router(idx1:end-idx2);
    [canal] = canalSurface(xx2,yy2,zz2,N, B, RRouter);
    
    values(1).Router = Router(idx1:end-idx2);
    values(1).xyz_distance = xyz_distance;
elseif strcmp(fit_type, 'ellipses')    
    RR1 = R1(idx1:end-idx2);
    RR2 = R2(idx1:end-idx2);
    alpha = alpha(idx1:end-idx2);
    canal = canalSurface_ellipse(xx2, yy2, zz2, N, B, RR1, RR2, alpha);
end

%Plot the canal surface
if plotSurface
    if ~isa(handles, 'struct') %plots in a normal figure
        plot_surface(xx2, yy2, zz2, data, canal, canalview, set_num, ...
            numDemos, fit_type, plot_method);
    else %plots in the GUI
        plot_surface2(xx2, yy2, zz2, data, canal, canalview, numDemos, ...
            plot_method, handles)
    end
end

%Populate remaining values into output cell
values(1).xmean = xx2;
values(1).ymean = yy2;
values(1).zmean = zz2;
values(1).N = N';
values(1).B = B';
values(1).T = T';
end

function plot_boundaries(data, B, xyz_distance, numDemos, threshold)
%Plots the boundaries in x, y, and z direction

%Plots the mean, trajectories, and boundary with an extra threshold
figure;
subplot(2,3,1);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,1),'k');
end
plot(B(2).xmean, 'b');
plot(xyz_distance(:,1)+B(2).xmean + threshold, 'r');
plot(-xyz_distance(:,1)+B(2).xmean - threshold, 'r');
title(sprintf('Set %i - X', set_num));axis square;grid;
hold off;

subplot(2,3,2);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,2),'k');
end
plot(B(2).ymean, 'b');
plot(xyz_distance(:,2)+B(2).ymean + threshold, 'r');
plot(-xyz_distance(:,2)+B(2).ymean - threshold, 'r');
title(sprintf('Set %i - Y', set_num));axis square;grid;
hold off;

subplot(2,3,3);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,3),'k');
end
plot(B(2).zmean, 'b');
plot(xyz_distance(:,3)+B(2).zmean + threshold, 'r');
plot(-xyz_distance(:,3)+B(2).zmean - threshold, 'r');
title(sprintf('Set %i - Z', set_num));axis square;grid;
hold off;

%Plots the data - mean and boundary
subplot(2,3,4);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,1)-B(2).xmean,'k');
end
plot(xyz_distance(:,1) + threshold, 'b');
title('Bound X');axis square;grid;
hold off;

subplot(2,3,5);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,2)-B(2).xmean,'k');
end
plot(xyz_distance(:,2) + threshold, 'b');
title('Bound Y');axis square;grid;
hold off;

subplot(2,3,6);
hold on;
for ii=1:numDemos
    plot(data{ii}(:,3)-B(2).xmean,'k');
end
plot(xyz_distance(:,1) + threshold, 'b');
title('Bound Z');axis square;grid;
hold off;
end

function plot_surface(xx2, yy2, zz2, data, canal, canalview, ...
    set_num, numDemos, fit_type, plot_method)
%Plots the canal surface

figure; hold on;
if strcmp('circles', plot_method)
    for kk = 1:size(canal,3)
        C = canal(:,:,kk);
        plot3(C(1,:),C(2,:),C(3,:),'k');
    end   
else
    surf(squeeze(canal(1,:,:)), squeeze(canal(2,:,:)), ...
        squeeze(canal(3,:,:)));
    shading interp;
    alpha 0.5;
    colormap bone;
end

plot3(xx2,yy2,zz2,'b','linewidth',2);

for ii=1:numDemos
    plot3(data{ii}(:,1), data{ii}(:,2), data{ii}(:,3),'r',...
        'linewidth',2);
end
title(sprintf('Set %i - Canal Surface using %s', ...
    set_num, fit_type));
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal;
view(canalview);
set(gcf, 'Position', get(0, 'Screensize'));
hold off;
end

function plot_surface2(xx2, yy2, zz2, data, canal, canalview, numDemos, ...
    plot_method, handles)
%Plots the canal surface for GUI

axes(handles.axes1);
cla reset;
hold on;
if strcmp('circles', plot_method)
    for kk = 1:size(canal,3)
        C = canal(:,:,kk);
        plot3(C(1,:),C(2,:),C(3,:),'k');
    end
else
    surf(squeeze(canal(1,1:2:end,1:2:end)), squeeze(canal(2,1:2:end,1:2:end)), ...
        squeeze(canal(3,1:2:end,1:2:end)));
    shading interp;
    alpha 0.5;
    colormap bone;
end

plot3(xx2,yy2,zz2,'b','linewidth',2);

for ii=1:numDemos
    plot3(data{ii}(:,1), data{ii}(:,2), data{ii}(:,3),'r',...
        'linewidth',2);
end
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal;
view(canalview);
hold off;
end
