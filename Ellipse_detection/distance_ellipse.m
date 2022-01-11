function dist = distance_ellipse(paramsC_i,Pts)
% Compute ellipse matrice from its parameters
    C = param2ellipse(paramsC_i);
    % If the product is negative then the point is inside the hull
    dist=0;
    for i = 1:size(Pts,2)
        if ((Pts(:,i)'*C*Pts(:,i)) < 0)
            dist = dist +1;
        end
    end
end