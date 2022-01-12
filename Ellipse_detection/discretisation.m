%% Parametre
R =10;
nb_points_disque = 30;

%% discretisation
nx =100;
ny =100;
nTheta = 100;
x = linspace(-sqrt(2)*R,sqrt(2)*R,nx);
y = linspace(-sqrt(2)*R,sqrt(2)*R,ny);
theta = linspace(0,2*pi,ntheta);

Res = zeros(100,100,100);
for i =1:100
    for j = 1:100
         for k = 1:100
                paramC1 = [0,0,R,sqrt(2)*R,0];
                paramC2 = [x(i),y(j),R,sqrt(2)*R,theta(k)];
                Pts1 = ellipsepoint(paramC1,nb_points_disque);
                Pts2 = ellipsepoint(paramC2,nb_points_disque);
                Res(i,j,k) = max(distance_ellipse(paramC1,Pts2),distance_ellipse(paramC2,Pts1));
         end
    end
    i
end

save('discretisation.mat');