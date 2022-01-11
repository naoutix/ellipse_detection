clear all
close all

display_ = 1;

% ellipse parameters xc, yc : center coordinates)
%                    Rx, Ry : semi-axis length
%                    theta : angle

% Create a random ellipse whose xc and theta are equal to 0
paramsC = [rand,rand,rand,rand,rand];

% Create random points in [-1,1]^2
nPoints = 200;
point = [-1+2*rand(1,nPoints);-1+2*rand(1,nPoints);ones(1,nPoints)];

if display_
    figure;
    plotellipse(paramsC); hold on;
    axis equal
end

% Compute ellipse matrice from its parameters
C = param2ellipse(paramsC);

% If the product is negative then the point is inside the hull
for i = 1:nPoints
    isInside(i) = (point(:,i)'*C*point(:,i)) < 0;
end

if display_
    % Collect points inside the elliptical hull
    for i=1:nPoints
        if isInside(i)
            plot(point(1,i),point(2,i),'c+');
        else
            plot(point(1,i),point(2,i),'r+');
        end
    end
end
