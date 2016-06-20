function [canal1,canal2,T,N,B,N2,B2] = canalSurface2(t,x,y,z,R,...
    vect_a,vect_b,plotIt,filterIt, cross_section_method)
% -----------------------------------------------------------------------
% A function to calculate a canal surface and its corresponding TNB
%
% Inputs:
%   t: time or interval
%   x,y,z: directrix or backbone curve of the canal
%   R: radii or radius function
%   vect_a, vect_b: set of orthonormal vectors for the cross-section
%   plotIt: true for plotting the canal and TNB false otherwise
%   filterIt: true for filtering the TNB false otherwise
%   cross_section_method: 1 or 2 depending on whether TNB frames are used
%                         to calculate cross-section alignment or the
%                         cross-section itself
%
% Output:
%   T,N,B: TNB vectors calculated using the first method
%   T,N2,B2: TNB vectors calculated using the second method
%   canal1: canal surface calculated using T,N,B
%   canal2: canal surface calculated using T,N2,B2
%
% -----------------------------------------------------------------------
% Code: Reza Ahmadzadeh 2016 (reza.ahmadzadeh@gatech.edu)
% April-6-2016
% Modified by: Roshni kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------
%% assessing data
% ----- check to see if the vector needs to be transposed -----
if size(t,1) > 1
    t = t.';
    x = x.';
    y = y.';
    z = z.';
    R = R.';
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
    if length(polyder(cx)) < 3; cxd(ii,2:end) = polyder(cx);
    else cxd(ii,:) = polyder(cx);
    end
    cxdd(ii,:) = polyder(cxd(ii,:));
end
% -----> in y direction
cyd = zeros(spy.pieces,spy.order-1);
cydd = zeros(spy.pieces,spy.order-2);
for ii = 1:spy.pieces
    cy = spy.coefs(ii,:);
    if length(polyder(cy)) < 3; cyd(ii,2:end) = polyder(cy);
    else cyd(ii,:) = polyder(cy);
    end
    cydd(ii,:) = polyder(cyd(ii,:));
end
% -----> in z direction
czd = zeros(spz.pieces,spz.order-1);
czdd = zeros(spz.pieces,spz.order-2);
for ii = 1:spz.pieces
    cz = spz.coefs(ii,:);
    if length(polyder(cz)) < 3; czd(ii,2:end) = polyder(cz);
    else czd(ii,:) = polyder(cz);
    end
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
N = zeros(size(rp));
for ii = 1:size(rp,2)
    g = norm(rp(:,ii),2);
    T(:,ii) = rp(:,ii) / g;
    gp = dot(rp(:,ii),rpp(:,ii))/g;
    kappa = norm(cross(rpp(:,ii),rp(:,ii)),2)/g^3;
    n = (g*rpp(:,ii)-gp*rp(:,ii))/(g^3*kappa);
    N(:,ii) = n / norm(n,2);
end
B = cross(T,N);
for ii = 1:size(B,2)
    B(:,ii) = B(:,ii) / norm(B(:,ii),2);
end

if filterIt
    a = 1;
    b = ones(1,20)/20;
    T = filter(b,a,T);
    B = filter(b,a,B);
    N = filter(b,a,N);
end

% ===== Second way of calculating TNB =====
% this method gives less noisier results (but the N and B vectors are
% rotated in their plane)
% TODO: the random vector can be a vector new to the first vector of T
randomVector = rand(3,1);
randomVector = randomVector / norm(randomVector,2);
N2 = cross(T,repmat(randomVector,1,size(T,2)));
% if we want to filter it has to be done before the normalization step
if filterIt
    a = 1;
    b = ones(1,20)/20;
    N2 = filter(b,a,N2);
end
for ii = 1:size(N2,2)
    N2(:,ii) = N2(:,ii) / norm(N2(:,ii),2);
end
B2 = cross(T,N2);
if filterIt
    B2 = filter(b,a,B2);
end
for ii = 1:size(B2,2)
    B2(:,ii) = B2(:,ii) / norm(B2(:,ii),2);
end
%% plotting the results
if plotIt
    stepc = 10;
    figure;
    subplot(2,2,1);hold on
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),...
        T(1,1:stepc:end),T(2,1:stepc:end),T(3,1:stepc:end));
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),...
        N(1,1:stepc:end),N(2,1:stepc:end),N(3,1:stepc:end));
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),...
        B(1,1:stepc:end),B(2,1:stepc:end),B(3,1:stepc:end));
    axis equal
    subplot(2,2,2);hold on
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),...
        T(1,1:stepc:end),T(2,1:stepc:end),T(3,1:stepc:end));
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),...
        N2(1,1:stepc:end),N2(2,1:stepc:end),N2(3,1:stepc:end));
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),...
        B2(1,1:stepc:end),B2(2,1:stepc:end),B2(3,1:stepc:end));
    axis equal
    subplot(2,2,3);hold on
end

tc = 0:0.01:1; % arc-length on the circle
stepc = fix(0.01*length(x)); % 10% of the data is used
L = length(1:stepc:xn);

allC1 = zeros(3,length(tc),L);
k = 1;
for ii = 1:stepc:length(xn)
    C = repmat([xn(ii);yn(ii);zn(ii)],1,length(tc)) + R(1,ii)*...
        (N(:,ii) * cos(2*pi*tc) + B(:,ii) * sin(2*pi*tc));
    allC1(:,:,k) = C;
    k = k + 1;
    if plotIt
        plot3(C(1,:),C(2,:),C(3,:),'k');
    end
end
if plotIt
    plot3(x,y,z,'linewidth',2);
    axis equal
    subplot(2,2,4);hold on
end

allC2 = zeros(3,length(tc),L);
k = 1;
%Align cross-sections to the orthonormal vectors
if cross_section_method == 2
    for ii = 1:stepc:length(xn)
        C = repmat([xn(ii);yn(ii);zn(ii)],1,length(tc)) + ...
            R(1,ii)*(vect_a(ii,:)'*cos(2*pi*tc)...
            +vect_b(ii,:)'*sin(2*pi*tc));
        allC2(:,:,k) = C;
        k = k + 1;
        if plotIt
            plot3(C(1,:),C(2,:),C(3,:),'k');
        end
    end
%Align cross-sections to the TNB frame
else
    for ii = 1:stepc:length(xn)
        C = repmat([xn(ii);yn(ii);zn(ii)],1,length(tc))...
            +R(1,ii)*(N2(:,ii)* cos(2*pi*tc) + B2(:,ii) * sin(2*pi*tc));
        allC2(:,:,k) = C;
        k = k + 1;
        if plotIt
            plot3(C(1,:),C(2,:),C(3,:),'k');
        end
    end
end

if plotIt
    plot3(x,y,z,'linewidth',2);
    axis equal
end

canal1 = allC1;
canal2 = allC2;
clear allC1 allC2 stepc k C
end