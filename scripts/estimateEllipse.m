function [a,b,aV,bV,alpha] = estimateEllipse(P, C)
% input:
%       P: a set of points in a plane size = (N,2)
%
% output:
%       a,b: major and minor radius
%       aV,bV: major and minor unit vectors 
%       alpha: rotation of the ellipse w.r.t x axis
%
% this function uses two different approaches to enclose the points inside
% the smallest feasible ellipse.

if size(P, 1) ~= 2
    P = P';
end
if size(C, 1) ~= 2
    C = C';
end

tol = 1e-3;
r = 1;

N = size(P,2);
D = sqrt(sum((P - repmat(C,1,N)).^2));  % all the distances from center
[Ds,Idx] = sort(D,'ascend');            % sort the distances
aV = [(P(1,Idx(end))-C(1,1)) ; (P(2,Idx(end))-C(2,1))];
aV = aV / norm(aV,2);                   % normal vector to the major axis
c = aV(1,1)*C(1,1) + aV(2,1)*C(2,1);    % get the line perpendicular to eV
xl = linspace(-1,1,20);
yl = (c - aV(1,1)*xl)/aV(2,1);
bV = [xl(1,1)-C(1,1); yl(1,1)-C(2,1)]; 
bV = bV/norm(bV,2);                   % unit vector of the minor axis
alpha = atan(aV(2,1)/aV(1,1));          % get angle w.r.t x axis
a = Ds(end);                            % major axis of the ellipse
for ii = 1:N
    b1 = Ds(ii);                       % consider it as minor axis
    % test to see if the current ellipse covers the whole points
    isInEllipse = ((cos(alpha)*(P(1,:) - repmat(C(1,1),1,N)) + sin(alpha)*(P(2,:)-repmat(C(2,1),1,N))).^2)/(a^2) + ...
        ((sin(alpha)*(P(1,:) - repmat(C(1,1),1,N)) - cos(alpha)*(P(2,:)-repmat(C(2,1),1,N))).^2)/(b1^2);
%     cond = sum(round((isInEllipse <= 1)*r)/r);
    cond = sum(isInEllipse-1 < tol);
    if (cond < N+tol && cond > N-tol)
        break;
    end
end

S = linspace(Ds(1),Ds(end),N*10);
for ii = 1:length(S)
    b2 = S(ii);
    % test to see if the current ellipse covers the whole points
    isInEllipse = ((cos(alpha)*(P(1,:) - repmat(C(1,1),1,N)) + sin(alpha)*(P(2,:)-repmat(C(2,1),1,N))).^2)/(a^2) + ...
        ((sin(alpha)*(P(1,:) - repmat(C(1,1),1,N)) - cos(alpha)*(P(2,:)-repmat(C(2,1),1,N))).^2)/(b2^2);
%     cond = sum(round((isInEllipse <= 1)*r)/r);
    cond = sum(isInEllipse-1 < tol);
    if (cond < N+tol && cond > N-tol)
        break;
    end
end

A1 = pi*a*b1;
A2 = pi*a*b2;
if A1 < A2
    b = b1;
else
    b = b2;
end

end