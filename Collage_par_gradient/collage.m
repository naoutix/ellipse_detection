function u = collage(r,s,interieur)
%% Calcul de l'imagette résultat
u=r;
r = double(r);
s = double(s);

%% Constante
[nb_lignes,nb_colonnes,nb_canaux] = size(s);
nb_pixels = nb_lignes * nb_colonnes;

%% Operateur gradient
e = ones(nb_pixels,1);
Dx = spdiags([-e e],[0 nb_lignes],nb_pixels,nb_pixels);
Dx(end-nb_lignes+1:end,:) = 0;
Dy = spdiags([-e e],[0 1],nb_pixels,nb_pixels);
Dy(nb_lignes:nb_lignes:end,:) = 0;

%% Matrice A 
A = -Dx'*Dx - Dy'*Dy;
bord_r = ones(nb_lignes,nb_colonnes);
bord_r(2:(nb_lignes-1),2:(nb_colonnes-1))=zeros(nb_lignes-2,nb_colonnes-2);

n_bord_r = (2*nb_lignes+2*nb_colonnes)-4;
indices_bord_r = find( bord_r == 1);

A(indices_bord_r,:) = sparse(1:n_bord_r,indices_bord_r,ones(n_bord_r,1),n_bord_r,nb_pixels);
for rgb=1:nb_canaux
    %% Mise en forme des donnees
    r_k= r(:,:,rgb);
    s_k= s(:,:,rgb);
    %% Calcul des gradients
    [gx,gy] =gradient(r_k);
    [grad_s_x,grad_s_y] = gradient(s_k);

    %% G
    gx(interieur)= grad_s_x(interieur);
    gy(interieur)= grad_s_y(interieur);
    gx = reshape(gx,[nb_pixels 1]);
    gy = reshape(gy,[nb_pixels 1]);   
    %% B
    b = Dx*gx+Dy*gy;
    b(indices_bord_r) = r_k(indices_bord_r);
    
    %% Resolution
    u_k = A\b;
    
    u(:,:,rgb)= reshape(u_k,[nb_lignes nb_colonnes]);
end
