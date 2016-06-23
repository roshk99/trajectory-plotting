function [T, N, B] = calculateTNB(t,xyz_mean)
% -----------------------------------------------------------------------
% A function that takes calculates a TNB frame from a mean trajectory
%
% Inputs:
%   t: time vector (0:end-1 of the data)
%   xyz_mean: contains the x,y,z components of the trajectory in the column
%
% Output:
%   T, N, B: contain unit vectors as columns for each time step (TNB
%   frames)
%
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------

x = xyz_mean(:,1)'; y = xyz_mean(:,2)'; z = xyz_mean(:,3)';

if size(t,1) > 1
    t = t.';
end

if size(x,1) > 1
    x = x.';
    y = y.';
    z = z.';
end

% ----- fit a spline to each dimension of the data -----
spx = spline(t,x);
spy = spline(t,y);
spz = spline(t,z);

% ----- evaluate the splines to make sure everything is alright -----
xn = ppval(spx,t);
yn = ppval(spy,t);
zn = ppval(spz,t);

% ------ get the coefficient for the derivative polynomials -----
% -----> in x direction
cxd = zeros(spx.pieces,spx.order-1);
cxdd = zeros(spx.pieces,spx.order-2);
for ii = 1:spx.pieces
    cx = spx.coefs(ii,:);
    if length(polyder(cx)) < 3; cxd(ii,2:end) = polyder(cx); else cxd(ii,:) = polyder(cx); end
    cxdd(ii,:) = polyder(cxd(ii,:));
end
% -----> in y direction
cyd = zeros(spy.pieces,spy.order-1);
cydd = zeros(spy.pieces,spy.order-2);
for ii = 1:spy.pieces
    cy = spy.coefs(ii,:);
    if length(polyder(cy)) < 3; cyd(ii,2:end) = polyder(cy); else cyd(ii,:) = polyder(cy); end
    cydd(ii,:) = polyder(cyd(ii,:));
end
% -----> in z direction
czd = zeros(spz.pieces,spz.order-1);
czdd = zeros(spz.pieces,spz.order-2);
for ii = 1:spz.pieces
    cz = spz.coefs(ii,:);
    if length(polyder(cz)) < 3; czd(ii,2:end) = polyder(cz); else czd(ii,:) = polyder(cz); end
    czdd(ii,:) = polyder(czd(ii,:));
end

% ----- make the polynomials using the coefficients -----
% -----> First and second derivatives of the spline.
spxd = mkpp(spx.breaks,cxd);
spyd = mkpp(spy.breaks,cyd);
spzd = mkpp(spz.breaks,czd);

spxdd = mkpp(spx.breaks,cxdd);
spydd = mkpp(spy.breaks,cydd);
spzdd = mkpp(spz.breaks,czdd);

% ----- Evaluate and get the real numbers for the derivatives -----
xd = ppval(spxd,t);
yd = ppval(spyd,t);
zd = ppval(spzd,t);

xdd = ppval(spxdd,t);
ydd = ppval(spydd,t);
zdd = ppval(spzdd,t);

%% calculate TNB
rp = [xd;yd;zd];
rpp = [xdd;ydd;zdd];
% ===== first way of calculating TNB =====
T = zeros(size(rp));
for ii = 1:size(rp,2)
    g = norm(rp(:,ii),2);
    T(:,ii) = rp(:,ii) / g;
end

randomVector = rand(3,1); % TODO: the random vector can be a vector new to the first vector of T
randomVector = randomVector / norm(randomVector,2);
N = cross(T,repmat(randomVector,1,size(T,2)));
% if we want to filter it has to be done before the normalization step
for ii = 1:size(N,2)
    N(:,ii) = N(:,ii) / norm(N(:,ii),2);
end
B = cross(T,N);
for ii = 1:size(B,2)
    B(:,ii) = B(:,ii) / norm(B(:,ii),2);
end
end
