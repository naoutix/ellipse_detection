close all
home
domain_bounds=[50,50];
data_points=[20 20;20 30;30 20;30 30];
%another possible example
%punto_dato=[25 25];
use_plot=1;
T=fast_marching(data_points,domain_bounds,use_plot);
mesh(T);