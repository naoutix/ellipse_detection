clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Paramètres :
beta = 1.0;                                     % Facteur de distance
S =130;                                         % seuil de bon postionnement 
gamma = 5.0;                                    % parametre sigmoide
T =0.1;                                         % Temperature
lambda = 100.0;                                 % facteur de poisson
alpha =0.99;                                    % coefficient de croissance
rose = [253 108 158]/255;
k_max = 5000;                                   % nb d'iteration max

N=10;                                            % Nb d'ellipse initaux

nb_points_disque = 30;                          % Affichage elipse
increment_angulaire = 2*pi/nb_points_disque;    % increment d'affichage
theta = 0:increment_angulaire:2*pi;             % Liste des thetas


nb_affichages = 1000;                               % affichage frame
pas_entre_affichages = floor(k_max/nb_affichages);
temps_pause = 0.0005;
%% Lecture et affichage de l'image :
I = imread('colonie.png');
I = rgb2gray(I);
I = double(I);
I = I(1:400,100:450);
figure('Name',['Detection a partir ' num2str(N) ' flamants roses'],'Position',[0.25*L,0,0.75*L,0.5*H]);
[nb_lignes,nb_colonnes] = size(I);

%% Tirage aléatoire d'une configuration initiale et calcul des niveaux de gris moyens :
c = zeros(N,5);
I_moyen_disques = zeros(N,1);
U_disques = zeros(N,1);
for i = 1:N
    pos = nb_colonnes*rand(1,2);
    a= rand*20+5;
    b= betarnd(4,4)*2*a;
    angle = rand*2*pi;
    c_i = [pos,rand*2*pi,a,b];
	c(i,:) = c_i;
	I_moyen_disques(i) = calcul_I_moyen_ellipse_variable(I,c_i);
    somme=0;
    for j = 1:N
        dist = sqrt((c(i,1)-c(j,1))^2 + (c(i,2)-c(j,2))^2);
        if dist <= sqrt(2)*(c(i,5)+c(j,5))/2
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
			x_prime = c(j,4)*cos(theta);
			y_prime = c(j,5)*sin(theta);
            angle = c(j,3);
            x_affich = x_prime*cos(angle)+y_prime*sin(angle);
            y_affich = -x_prime*sin(angle)+y_prime*cos(angle);
            x_affich = c(j,1)+x_affich;
            y_affich = c(j,2)+y_affich;
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
    
    for i=1:nouvN
        % Tirage aléatoire d'un nouveau disque et calcul du niveau de gris moyen :
        a= rand*20+5;
        b= betarnd(4,4)*2*a;
        c_i = [nb_colonnes*rand nb_lignes*rand,rand*2*pi,a,b];
        c = [c; c_i];
        I_moyen_disques = [I_moyen_disques ; calcul_I_moyen_ellipse_variable(I,c_i)];
    end
    
 
    
    %% Calcul de l'energie ind
    U_i = zeros(size(c,1),1);
    
    for i = 1:size(c,1)
        U_i(i) = 1 - (2/(1+exp(-gamma*(I_moyen_disques(i)/S-1))));
    end
    %% tri des energies
    [U_i, ind] = sort(U_i, 'desc');
    
    c = c(ind,:);
    I_moyen_disques = I_moyen_disques(ind,:);
    
    %% morts des cercles
    i=1;
    while i <= size(c,1)
        somme =0;
        sum_Ci=0;
       
        for j=1:size(c,1)
            for l=1:size(c,1)
                if j ~= l
                    dist = sqrt((c(j,1)-c(l,1))^2 + (c(j,2)-c(l,2))^2);
                    if dist<= sqrt(2)*(c(i,5)+c(j,5))/2
                        somme = somme +1;
                        if j~=i && l~=i
                            sum_Ci = sum_Ci +1;
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
			x_prime = c(j,4)*cos(theta);
			y_prime = c(j,5)*sin(theta);
            angle = c(j,3);
            x_affich = x_prime*cos(angle)+y_prime*sin(angle);
            y_affich = -x_prime*sin(angle)+y_prime*cos(angle);
            x_affich = c(j,1)+x_affich;
            y_affich = c(j,2)+y_affich;
        indices = find(x_affich>0 & x_affich<nb_colonnes & y_affich>0 & y_affich<nb_lignes);
        plot(x_affich(indices),y_affich(indices),'Color',rose,'LineWidth',3);
        title(['Detection en cours : ' num2str(size(c,1)) ' flamants roses trouvés']);
    end
    pause(temps_pause);
	

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
end