clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

load('discretisation.mat')
%% Paramètres :
R = 10;                                         % Rayon des disques 10 ( dans le cas ellipitique a = R ,b = sqrt(2)*R)
beta = 1.0;                                     % Facteur de prise en compte des distances
S =130;                                         % seuil de bon positionnement 
gamma = 5.0;                                    % parametre sigmoide
T =0.1;                                         % Temperature
lambda = 50.0;                                 % facteur de poisson
alpha =0.99;                                    % coefficient de decroissance
N=10;                                            % Nb de cercle initaux
k_max = 5000;                                   % nb d'iteration max

nb_points_disque = 30;                          % precision des ellipse
Pourcentage_point_interne_tol = 25;             % Pourcentage de chevauchement des ellipses
rose = [253 108 158]/255;                       % couleur 


nb_affichages = 1000;                               % affichage frame
pas_entre_affichages = floor(k_max/nb_affichages);
temps_pause = 0.0005;

distmax =  round(nb_points_disque * Pourcentage_point_interne_tol/100);  % Nombre de point maximal de chevauchement
%% Lecture et affichage de l'image :
I = imread('colonie.png');
I = rgb2gray(I);
I = double(I);
I = I(1:400,100:450);
%I = imbinarize(I);
figure('Name',['Detection de ' num2str(N) ' flamants roses'],'Position',[0.25*L,0,0.75*L,0.5*H]);
[nb_lignes,nb_colonnes] = size(I);

%% Tirage aléatoire d'une configuration initiale et calcul des niveaux de gris moyens :
c = zeros(N,3);
Points = zeros(3,nb_points_disque,N);
I_moyen_disques = zeros(N,1);
U_disques = zeros(N,1);
for i = 1:N
    pos = nb_colonnes*rand(1,2);
    c_i = [pos,rand*2*pi];
	c(i,:) = c_i;
	I_moyen_disques(i) = calcul_I_moyen_ellipse(I,c_i,R);
    somme=0;
    paramsC_i = [c(i,1),c(i,2),R,sqrt(2)*R,c(i,3)];
    PtsC_i = ellipsepoint(paramsC_i,nb_points_disque);
    Points(:,:,i)=PtsC_i;
    for j = 1:N
        paramsC_j = [c(j,1),c(j,2),R,sqrt(2)*R,c(j,3)];
        %dist = max(distance_ellipse(paramsC_j,PtsC_i),distance_ellipse(paramsC_i,ellipsepoint(paramsC_j,nb_points_disque)));
        x_interpol = round((paramsC_i(1)-paramC_j(1)-x(1))/((x(nx)-x(1))/nx)) +1;
        y_interpol = round((paramsC_i(2)-paramC_j(2)-y(1))/((y(ny)-y(1))/ny)) +1;
        theta_interpol = round((abs(paramsC_i(5)-paramC_j(5))-theta(1))/((theta(nTheta)-theta(1))/nTheta)) +1;
        dist = Res(x_interpol,y_interpol,theta_interpol);
        if dist > distmax
           somme = somme + 1;
        end
    end
    U_disques(i) = beta*somme - I_moyen_disques(i);
end
liste_k = 0;
liste_U = 0;
I_moyen_config = mean(I_moyen_disques);
liste_I_moyen_config = I_moyen_config;

%% Affichage de la configuration initiale :
subplot(1,2,1);
imagesc(I);
axis image;
axis off;
colormap gray;
hold on;
for j = 1:N
            x_affich = Points(1,:,j);
            y_affich = Points(2,:,j);
			indices = find(x_affich>0 & x_affich<nb_colonnes & y_affich>0 & y_affich<nb_lignes);
			plot(x_affich(indices),y_affich(indices),'Color',rose,'LineWidth',3);
end
pause(temps_pause);

%% Courbe d'évolution du niveau de gris moyen :
subplot(1,2,2);
plot(liste_k,liste_U,'.','Color',rose);
axis([0 k_max/1000 -400 0]);
set(gca,'FontSize',20);
xlabel('Nombre d''iterations','FontSize',30);
ylabel('Energie','FontSize',30);

        
%% Recherche de la configuration optimale :
continuer = true;
k=1;
while continuer
    %% Tirage des nouveaux disques
    c_anc = c;
    nouvN = poissrnd(lambda);
    [nl,nc,N_anc] = size(Points);
    Points2 = zeros(nl,nc,N_anc+nouvN);
    Points2(:,:,1:N_anc)=Points;
    for i=1:nouvN
        % Tirage aléatoire d'un nouveau disque et calcul du niveau de gris moyen :
        c_i = [nb_colonnes*rand nb_lignes*rand,rand*2*pi];
        c = [c; c_i];
        paramsC_i = [c_i(1),c_i(2),R,sqrt(2)*R,c_i(3)];
        Points2(:,:,N_anc+i)=ellipsepoint(paramsC_i,nb_points_disque);
        I_moyen_disques = [I_moyen_disques ; calcul_I_moyen_ellipse(I,c_i,R)];
    end
    Points = Points2;
 
    
    %% Calcul de l'energie ind
    U_i = zeros(size(c,1),1);
    
    for i = 1:size(c,1)
        U_i(i) = 1 - (2/(1+exp(-gamma*(I_moyen_disques(i)/S-1))));
    end
    %% tri des energies
    [U_i, ind] = sort(U_i, 'desc');
    
    c = c(ind,:);
    Points = Points(:,:,ind);
    I_moyen_disques = I_moyen_disques(ind,:);
    
    %% morts des cercles
    i=1;
    while i <= size(c,1)
        somme =0;
        sum_Ci=0;
       
        for j=1:size(c,1)
            paramsC_j = [c(j,1),c(j,2),R,sqrt(2)*R,c(j,3)];
            for l=1:size(c,1)
                paramsC_l = [c(l,1),c(l,2),R,sqrt(2)*R,c(l,3)];
                if j ~= l
                    dist1 = sqrt((c(j,1)-c(l,1))^2 + (c(j,2)-c(l,2))^2);
                    if dist1 <= 2*sqrt(2)*R
                        %dist = max(distance_ellipse(paramsC_j,Points(:,:,l)),distance_ellipse(paramsC_l,Points(:,:,j)));
                            x_interpol = round((paramsC_i(1)-paramC_j(1)-x(1))/((x(nx)-x(1))/nx)) +1;
                            y_interpol = round((paramsC_i(2)-paramC_j(2)-y(1))/((y(ny)-y(1))/ny)) +1;
                            theta_interpol = round((abs(paramsC_i(5)-paramC_j(5))-theta(1))/((theta(nTheta)-theta(1))/nTheta)) +1;
                            dist = Res(x_interpol,y_interpol,theta_interpol);
                        if dist> distmax
                            somme = somme +1;
                            if j~=i && l~=i
                                sum_Ci = sum_Ci +1;
                            end
                        end
                    end 
                end
            end
        end  
        U_sans_i = sum(U_i) - U_i(i) + beta * sum_Ci;
        U = sum(U_i) + beta*somme;
        
        proba = (lambda) / (lambda* + exp((U_sans_i - U)/(T*alpha.^k)));
        if rand < proba
            c(i,:) =[];
            U_i(i) =[];
            Points(:,:,i)=[];
            I_moyen_disques(i) = [];
            i = i - 1;
        end
        i = i+1;
    end
   
    %% Verification du changement des cerles
    if ~isequal(c_anc,c)
        T=alpha*T;
        lambda=alpha*lambda;
    else
        continuer = false;
    end
    
    %% Affichage cercles
    hold off;
    subplot(1,2,1);
    imagesc(I);
    axis image;
    axis off;
    colormap gray;
    hold on;
    for j = 1:size(c,1)
            x_affich = Points(1,:,j);
            y_affich = Points(2,:,j);
        indices = find(x_affich>0 & x_affich<nb_colonnes & y_affich>0 & y_affich<nb_lignes);
        plot(x_affich(indices),y_affich(indices),'Color',rose,'LineWidth',3);
        title(['Detection en cours : ' num2str(size(c,1)) ' flamants roses trouvés']);
    end
	

	% Courbe d'évolution du niveau de gris moyen :
	if rem(k,pas_entre_affichages)==0
		liste_k = [liste_k k];
		I_moyen_config = mean(I_moyen_disques);
		liste_I_moyen_config = [liste_I_moyen_config I_moyen_config];
        liste_U = [liste_U U];
		subplot(1,2,2);
		plot(liste_k,liste_U,'.-','Color',rose,'LineWidth',3);
		axis([0 max(k_max/1000,1.05*k) -400 0]);
		set(gca,'FontSize',20);
		xlabel('Nombre d''iterations','FontSize',30);
		ylabel('Energie','FontSize',30);
    end
    k=k+1;
    disp("Nombre de flamants")
    disp(size(c,1))
    pause(temps_pause);
end