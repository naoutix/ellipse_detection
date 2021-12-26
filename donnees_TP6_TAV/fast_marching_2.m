function D=fast_marching_2(frontiere,size)
% 0 - calculé
% 1 - proche
% 2 - restant

nb_point_frontiere = length(frontiere(:,1));
voisin=[1 0;-1 0;0 1;0 -1];
% Distance infinie
D=ones(size(1),size(2))*inf;
labels=2*ones(size(1),size(2));

%Mise en place des points calcule
lin_ind=sub2ind(size,frontiere(:,1),frontiere(:,2));
D(lin_ind)=0;
labels(lin_ind)=0;

%% intialisation
for i=1:nb_point_frontiere
    for k=[0,1,-1]
        for l=[0,1,-1]
            x_k =frontiere(i,1)+k;
            y_k =frontiere(i,2)+l;
            % test des limite
            if x_k > 0 && x_k <=size(1) && y_k >0 && y_k<=size(2)
                % test deja calculé
                p_voisin = sub2ind(size,x_k,y_k);
                if labels(p_voisin) ~= 0
                    labels(p_voisin) = 1;
                    %D(p_voisin)=sqrt(k^2+l^2);
                end
            end
        end
    end
end
 
%% fast marching
stop=0;
while ~stop
    proche=find(labels==1);
    if ~isempty(proche)
        [~,index]=sort(D(proche));
        proche_min = proche(index);
        for pixel_proche=1:length(proche)
            p1 = proche_min(pixel_proche);

            %passage a calculé
            labels(p1)=0;

            [coord_x,coord_y]=ind2sub(size,p1);
            pointP1=[coord_x,coord_y];
            min_dist = inf;
            for k=[0,1,-1]
                for l=[0,1,-1]
                    x_k =pointP1(1)+k;
                    y_k =pointP1(2)+l;
                    % test des limite
                    if x_k > 0 && x_k <=size(1) && y_k >0 && y_k<=size(2)
                        % test deja calculé
                        p_voisin = sub2ind(size,x_k,y_k);
                        if labels(p_voisin) == 0
                            d_voisin = sqrt(k^2+l^2) + D(p_voisin);
                            if d_voisin < min_dist
                                min_dist = d_voisin;
                            end
                        elseif labels(p_voisin) == 2
                                labels(p_voisin)=1;
                        end
                    end
                end 
             end
            D(p1) = min_dist;
        end
    else
        stop=1;
    end
end
