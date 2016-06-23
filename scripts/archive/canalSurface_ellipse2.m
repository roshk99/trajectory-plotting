function [canal] = canalSurface_ellipse2(xn, yn, zn, N2, B2 ,R1,...
    R2,alpha)

tc = 0:0.01:1; % arc-length on the circle
stepc = fix(0.01*length(xn)); % 10% of the data is used
L = length(1:stepc:xn);
allC1 = zeros(3,length(tc),L);
k = 1;
for ii = 1:stepc:length(xn)
    a = R1(ii); b = R2(ii); alpha_val = alpha(ii);
    
    x = a*cos(2*pi*tc)*cos(alpha_val) - b*sin(2*pi*tc)*sin(alpha_val);
    y = a*cos(2*pi*tc)*sin(alpha_val) + b*sin(2*pi*tc)*cos(alpha_val);
    C = repmat([xn(ii);yn(ii);zn(ii)],1,length(tc)) + N2(:,ii)*x + B2(:,ii)*y;
    allC1(:,:,k) = C;
    k = k + 1;
end

canal = allC1;
clear allC1 stepc k C
end