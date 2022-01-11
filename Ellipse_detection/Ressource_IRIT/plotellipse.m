function h = plotellipse(params,s,T)
% function [h, pts] = plotellipse(params,s,T)
%          params (required) = [uc, vc, Ru, Rv, theta]; (theta in radians)
%          s      (optional) = string for "line types, plot symbols and colors"
%          T      (optional) = homogeneous 3-by-3 transformation

switch nargin
          case 1
            s          =            '-';
            T          =       eye(3,3);
          case 2
            T          =       eye(3,3);
end    

t   = linspace(0,pi*2,1000);
x   = params(3) * cos(t);
y   = params(4) * sin(t);
nx  = x*cos(params(5))-y*sin(params(5)) + params(1); 
ny  = x*sin(params(5))+y*cos(params(5)) + params(2);

pts = T * augment([nx;ny]);
pts = pts(1:2,:);
h   = plot(pts(1,:),pts(2,:),s);
