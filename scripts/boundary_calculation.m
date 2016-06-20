function values = boundary_calculation(data, set_num, ...
    canalview, howmany, idx1, idx2, plotSurface, fit_type, canal_method)
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
%   canal_method: 1 or 2 depending on which method to create the canal is
%                 used
%
% Output:
%   values: cell with various canal surface parameters
%           values(1).xmean = mean in x direction (point_num x 1)
%           values(1).ymean = mean in y direction (point_num x 1)
%           values(1).zmean = mean in z direction (point_num x 1)
%           values(1).N2 = normal component of TNB frame (3 x point_num)
%           values(1).B2 = binormal component of TNB frame (3 x point_num)
%           values(1).T = tangent component of TNB frame (3 x point_num)
%
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------

%Initialize flags for internal plotting and filtering
plotBoundaries = false;
plotSurface_internal = false;
filterSurface_internal = false;

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

%Calculates the mean trajectory
B = calculate_mean(allXs, allYs, allZs);

%Finds the boundary values (radii and orientation vectors)
if strcmp(fit_type, 'circle')
    if canal_method == 1
        [Router, xyz_distance] = find_boundaries([B(2).xmean, ...
            B(2).ymean, B(2).zmean], allXs, allYs, allZs);
    else
        [Router, xyz_distance, vect_a, vect_b] = find_boundaries2([B(2).xmean, ...
            B(2).ymean, B(2).zmean], allXs, allYs, allZs);
    end
    
    if plotBoundaries
        plot_boundaries(data, B, xyz_distance, numDemos, threshold);
    end
elseif strcmp(fit_type, 'ellipse')
    [R1, R2, alpha, vect_a, vect_b] = find_boundaries_ellipse([B(2).xmean, ...
        B(2).ymean, B(2).zmean], allXs, allYs, allZs);
end


%Get only subset of values within specified boundary
t = linspace(0,1,nPoints).';
tt = t(idx1:end-idx2);
xx2 = B(2).xmean(idx1:end-idx2);
yy2 = B(2).ymean(idx1:end-idx2);
zz2 = B(2).zmean(idx1:end-idx2);

values = struct([]);

%Get the canal surface depending on the fit type
if strcmp(fit_type, 'circle')
    RRouter = Router(idx1:end-idx2);
    
    if canal_method == 1
        [~,canal,T,~,~,N2,B2] = ...
        canalSurface(tt,xx2,yy2,zz2,RRouter,plotSurface_internal,...
        filterSurface_internal);
    else
        vect_a = vect_a(idx1:end-idx2, :);
        vect_b = vect_b(idx1:end-idx2, :);
        [~,canal,T,~,~,N2,B2] = ...
            canalSurface2(tt,xx2,yy2,zz2,RRouter,vect_a, vect_b, ...
            plotSurface_internal,filterSurface_internal);
    end
    
    values(1).Router = Router(idx1:end-idx2);
    values(1).xyz_distance = xyz_distance;
elseif strcmp(fit_type, 'ellipse')
    RR1 = R1(idx1:end-idx2);
    RR2 = R2(idx1:end-idx2);
    Ralpha = alpha(idx1:end-idx2);
    Rvect_a = vect_a(idx1:end-idx2,:);
    Rvect_b = vect_b(idx1:end-idx2,:);
    [canal,T,~,~,N2,B2] = ...
        canalSurface_ellipse(tt,xx2,yy2,zz2,RR1,RR2,Ralpha,Rvect_a,Rvect_b,plotSurface_internal,...
        filterSurface_internal, canal_method);
end

%Plot the canal surface
if plotSurface
    plot_surface(xx2, yy2, zz2, data, canal, ...
        canalview, set_num, numDemos)
end

%Populate remainig values into output cell
values(1).xmean = xx2;
values(1).ymean = yy2;
values(1).zmean = zz2;
values(1).N2 = N2';
values(1).B2 = B2';
values(1).T = T';
end

function B = calculate_mean(allXs, allYs, allZs)
%Returns the mean trajectory
B = struct([]);
B(2).xmean = mean(allXs,2);
B(2).ymean = mean(allYs,2);
B(2).zmean = mean(allZs,2);
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
    set_num, numDemos)
%Plots the canal surface

figure;hold on;
for kk = 1:size(canal,3)
    C = canal(:,:,kk);
    plot3(C(1,:),C(2,:),C(3,:),'k');
end
plot3(xx2,yy2,zz2,'b','linewidth',2);
%     surf(squeeze(canal(1,:,:)), squeeze(canal(2,:,:)), squeeze(canal(3,:,:)));
%     shading interp;
%     alpha 0.5;
%     colormap bone;

for ii=1:numDemos
    plot3(data{ii}(:,1), data{ii}(:,2), data{ii}(:,3),'r',...
        'linewidth',2);
end
title(sprintf('Set %i - Canal Surface', set_num));
xlabel('X'); ylabel('Y'); zlabel('Z');
axis square
view(canalview);
set(gcf, 'Position', get(0, 'Screensize'));
hold off;
end
