function [canal,T,N,B,N2,B2] = canalSurface_ellipse(t,x,y,z,R1,...
    R2,alpha,vect_a,vect_b,plotIt,filterIt, canal_method)
%% assessing data
% ----- check to see if the vector needs to be transposed -----
if size(t,1) > 1
    t = t.';
    x = x.';
    y = y.';
    z = z.';
    R1 = R1.';
    R2 = R2.';
    alpha = alpha.';
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
randomVector = rand(3,1);           % TODO: the random vector can be a vector new to the first vector of T
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
    subplot(2,1,1);hold on
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),T(1,1:stepc:end),T(2,1:stepc:end),T(3,1:stepc:end));
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),N(1,1:stepc:end),N(2,1:stepc:end),N(3,1:stepc:end));
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),B(1,1:stepc:end),B(2,1:stepc:end),B(3,1:stepc:end));
    axis equal
    subplot(2,1,2);hold on
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),T(1,1:stepc:end),T(2,1:stepc:end),T(3,1:stepc:end));
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),N2(1,1:stepc:end),N2(2,1:stepc:end),N2(3,1:stepc:end));
    quiver3(x(1:stepc:end),y(1:stepc:end),z(1:stepc:end),B2(1,1:stepc:end),B2(2,1:stepc:end),B2(3,1:stepc:end));
    axis equal
end

tc = 0:0.01:1; % arc-length on the circle
stepc = fix(0.01*length(x)); % 10% of the data is used
L = length(1:stepc:xn);
alpha  = (alpha/pi*180 + 125)*pi/180;
allC1 = zeros(3,length(tc),L);
k = 1;
for ii = 1:stepc:length(xn)
    a = R1(ii); b = R2(ii); alpha_val = alpha(ii);
    
    x = a*cos(2*pi*tc)*cos(alpha_val) - b*sin(2*pi*tc)*sin(alpha_val);
    y = a*cos(2*pi*tc)*sin(alpha_val) + b*sin(2*pi*tc)*cos(alpha_val);
    if canal_method == 1
        C = repmat([xn(ii);yn(ii);zn(ii)],1,length(tc)) + N2(:,ii)*x + B2(:,ii)*y;
    else
        C = repmat([xn(ii);yn(ii);zn(ii)],1,length(tc)) + vect_a(ii,:)'*x + vect_b(ii,:)'*y;
    end
    allC1(:,:,k) = C;
    k = k + 1;
end

canal = allC1;
clear allC1 stepc k C
end