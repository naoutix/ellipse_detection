function resultat = calcul_I_moyen_ellipse(I,c_i,R)

[nb_lignes,nb_colonnes] = size(I);
abscisse = c_i(1);
ordonnee = c_i(2);
angle = c_i(3);
a = R;
b = sqrt(2)*R;
nb_pixels = 0;
somme_nvg = 0;
for j = max(1,floor(abscisse-b)):min(nb_colonnes,ceil(abscisse+b))
	for i = max(1,floor(ordonnee-b)):min(nb_lignes,ceil(ordonnee+b))
		abscisse_relative = j-abscisse;
		ordonnee_relative = i-ordonnee;
		if ((abscisse_relative*cos(angle)-ordonnee_relative*sin(angle))/a)^2+((abscisse_relative*sin(angle)+ordonnee_relative*cos(angle))/b)^2<=1
			nb_pixels = nb_pixels+1;
			somme_nvg = somme_nvg+I(i,j);
		end
	end
end
resultat = somme_nvg/nb_pixels;
