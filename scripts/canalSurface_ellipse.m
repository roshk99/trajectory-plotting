function [canal] = canalSurface_ellipse(xn, yn, zn, N, B ,R1,...
    R2,alpha)
% -----------------------------------------------------------------------
% A function that calculates the canal for elliptical cross-sections
%
% Inputs:
%   xn, yn, zn: the mean trajectory components
%   N, B: part of the TNB coordinate frames
%   R1, R2: the major and minor radii of the ellipse
%   alpha: the angle between the N-axis and and major axis of the ellipse
% Output:
%   canal: canal surface calculated with each circle as canal(:,:,i)
%
% -----------------------------------------------------------------------
% Code: Roshni Kaushik 2016 (roshni.s.kaushik@gmail.com)
% -----------------------------------------------------------------------

tc = 0:0.01:1; % arc-length on the circle
stepc = fix(0.01*length(xn)); % 10% of the data is used
L = length(1:stepc:xn);
allC1 = zeros(3,length(tc),L);
k = 1;
for ii = 1:stepc:length(xn)
    a = R1(ii); b = R2(ii); alpha_val = alpha(ii);
    
    x = a*cos(2*pi*tc)*cos(alpha_val) - b*sin(2*pi*tc)*sin(alpha_val);
    y = a*cos(2*pi*tc)*sin(alpha_val) + b*sin(2*pi*tc)*cos(alpha_val);
    C = repmat([xn(ii);yn(ii);zn(ii)],1,length(tc)) + N(:,ii)*x + B(:,ii)*y;
    allC1(:,:,k) = C;
    k = k + 1;
end

canal = allC1;
clear allC1 stepc k C
end