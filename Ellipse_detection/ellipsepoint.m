function Pts = ellipsepoint(paramsC_j,nb_points_disque)

    t   = linspace(0,pi*2,nb_points_disque);
    x   = paramsC_j(3) * cos(t);
    y   = paramsC_j(4) * sin(t);
    nx  = x*cos(paramsC_j(5))-y*sin(paramsC_j(5)) + paramsC_j(1); 
    ny  = x*sin(paramsC_j(5))+y*cos(paramsC_j(5)) + paramsC_j(2);
    Pts = [nx;ny;ones(1,length(nx))];
        
