function [canal] = canalSurface(xn, yn, zn, N, B, R)
% -----------------------------------------------------------------------
% A function that calculates the canal for circular cross-sections
%
% Inputs:
%   xn, yn, zn: the mean trajectory components
%   N, B: part of the TNB coordinate frames
%   R: the outer radius of the circle
%
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
    C = repmat([xn(ii);yn(ii);zn(ii)],1,length(tc))...
        +R(ii)*(N(:,ii)* cos(2*pi*tc) + B(:,ii) * sin(2*pi*tc));
    allC1(:,:,k) = C;
    k = k + 1;
end

canal = allC1;
clear allC1 stepc k C
end