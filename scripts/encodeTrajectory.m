function [values, canal] = encodeTrajectory(data, howmany, idx1, ...
    idx2, fit_type)
% -----------------------------------------------------------------------
% A function that takes smoothed data and generates a canal surface
%
% Inputs:
%   data: cell with each element point_numx3 vector of smoothed data
%         (output of parsetrajectory.m)
%   howmany: which trajectories to include in analysis (a list of indices)
%   idx1, idx2: values in the beginning and end of trajectories will be
%               eliminated when plotting the surface; data used will be
%               from idx1 to end-idx2
%   fit_type: 'circle', 'ellipse', or 'bspline' based on type of surface
%             desired
%
% Output:
%   values: cell with various canal surface parameters
%           values(1).xmean = mean in x direction (point_num x 1)
%           values(1).ymean = mean in y direction (point_num x 1)
%           values(1).zmean = mean in z direction (point_num x 1)
%           values(1).N = normal component of TNB frame (3 x point_num)
%           values(1).B = binormal component of TNB frame (3 x point_num)
%           values(1).T = tangent component of TNB frame (3 x point_num)
%   canal: stores the X,Y,Z values for plotting the canal surface
%
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------

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
    [Router, xyz_distance] = findBoundaries_circle(xyz_mean, allXs, allYs, ...
        allZs);
elseif strcmp(fit_type, 'ellipses')
    [R1, R2, alpha] = findBoundaries_ellipse(xyz_mean, ...
        allXs, allYs, allZs, T, N, B);
end


%Get only subset of values within specified boundary
xx2 = xyz_mean(idx1:end-idx2, 1);
yy2 = xyz_mean(idx1:end-idx2, 2);
zz2 = xyz_mean(idx1:end-idx2, 3);
N = N(:, idx1:end-idx2);
B = B(:, idx1:end-idx2);
values = struct([]);

%Get the canal surface depending on the fit type
if strcmp(fit_type, 'circles')
    RRouter = Router(idx1:end-idx2);
    canal = canalSurface_circle(xx2,yy2,zz2,N, B, RRouter);
    
    values(1).Router = RRouter;
    values(1).xyz_distance = xyz_distance(idx1:end-idx2,:);
elseif strcmp(fit_type, 'ellipses')    
    RR1 = R1(idx1:end-idx2);
    RR2 = R2(idx1:end-idx2);
    alpha = alpha(idx1:end-idx2);
    canal = canalSurface_ellipse(xx2, yy2, zz2, N, B, RR1, RR2, alpha);
end

%Populate remaining values into output cell
values(1).xmean = xx2;
values(1).ymean = yy2;
values(1).zmean = zz2;
values(1).N = N';
values(1).B = B';
values(1).T = T';
end

