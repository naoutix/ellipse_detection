% USAGE:	T=fast_marching(data_points);
%
% ARGUMENTS:
%				*data_points:	column vector with 2-dimensional input points.
%									i.e. data_points=[x1 y1; x2 y2; ...;xn yn].
%									All the coordinates should be integer values.
%				*siz:				size of the image domain. Domain is assumed to
%									be the rectangle [0 0 siz(1) siz(2)]. All the
%									input coordinates should be in that range.
%				*plot:			boolean (optional). enables or disables interactive plot.
% OUTPUT:
%				*T:			image with the distance map
function T=fast_marching(data_points,siz,plot)
%check arguments
% M1=max(data_points,1);
% M2=max(data_points,2);
% m1=min(data_points,1);
% m2=min(data_points,2);
if(sum(data_points(:,1)>siz(1)-1) || sum(data_points(:,1)<2) || sum(data_points(:,2)>siz(2)-1) || sum(data_points(:,2)<2))
    error('data points out of bounds');
end
% bug: estoy cambiando el valor de T en los alive
if(nargin<3)
    plot=1;
end
% tipos
% 0 - alive
% 1 - narrow band
% 2 - far away
dim=2;
delta=1;
largo_datos=size(data_points,1);
% hay que ver el tipo de conectividad
% estan ordenados en x,y,z,etc
vecinos=[1 0;-1 0;0 1;0 -1];
num_vecinos=size(vecinos,1);
% la mayoria infinito (los far away)
F=ones(siz(1),siz(2));
T=ones(siz(1),siz(2))*inf;
u=zeros(dim,1);
m=zeros(dim,1);
tipos=2*ones(siz(1),siz(2));
% pongo en cero los puntos inicialmente alive
lin_ind=sub2ind(siz,data_points(:,1),data_points(:,2));
T(lin_ind)=0;
tipos(lin_ind)=0;
% defino los narrow band
for k=1:largo_datos,
    for d=1:dim
        coord_vecinos(:,d)=data_points(k,d)+vecinos(:,d);
        ind=find( coord_vecinos(:,d)>siz(d) | coord_vecinos(:,d)<0 );
        coord_vecinos(ind,d)=data_points(k,d);
    end
    lin_ind=sub2ind(siz,coord_vecinos(:,1),coord_vecinos(:,2));
    tipos(lin_ind)=(tipos(lin_ind)~=0);
    T(lin_ind)=1./F(lin_ind);
end
% para cambiar el color del marker en un plot luego de hecho
% set(get(gca,'Children'),'Color',[1 0 0])
stop=0;
iter=1;
while ~stop,
	%fprintf('1');
	% busco cual de los NB tiene el T mas chico
	ind_narrow_band=find(tipos==1); % es un indice lineal
	if ~isempty(ind_narrow_band),
		[useless,ind_rel_min_T]=min(T(ind_narrow_band));
		ind_min_T=ind_narrow_band(ind_rel_min_T);
		% y lo paso a alive
		tipos(ind_min_T)=0;
		%  fprintf('2');
		% paso a los vecinos a la NB
		[coord_x,coord_y]=ind2sub(siz,ind_min_T);
		punto=[coord_x,coord_y];
		for d=1:dim
			coord_vecinos(:,d)=punto(d)+vecinos(:,d);
			ind=find(coord_vecinos(:,d)>siz(d) | coord_vecinos(:,d)<=0 );
			%REVISAR: hay que arreglar el borde
			coord_vecinos(ind,d)=punto(d)-(coord_vecinos(ind,d)-punto(d));
		end
		%fprintf('3');
		% busco los vecinos en far away y los paso a NB
		lin_ind=sub2ind(siz,coord_vecinos(:,1),coord_vecinos(:,2));
		ind_rel_far=find(tipos(lin_ind)==2);
		ind_far=lin_ind(ind_rel_far);
		tipos(ind_far)=1;
		% calculo el nuevo T, solo de los vecinos o de toda la narrow band???!!!
		% creo que es solo en los vecinos
		% primero calculo los m_i
		cant_vecinos_narrow=length(ind_far);
		%fprintf('4');
		%keyboard
		%for t=1:cant_vecinos_narrow,
		% hay que ver como hacer para recorrer solo los de la NB
		for t=1:4,
			current_vecino=coord_vecinos(t,:);
			F_curr=F(current_vecino(1),current_vecino(2));
			% calculo los vecinos del vecino (para actualizar su valor de T)
			%fprintf('5');
			for d=1:dim
				coord_vec_vec(:,d)=current_vecino(d)+vecinos(:,d);
				ind=find(coord_vec_vec(:,d)>siz(d) | coord_vec_vec(:,d)<=0 );
				% arreglo el borde (simetria, periodico)
				coord_vec_vec(ind,d)=current_vecino(d)-(coord_vec_vec(ind,d)-current_vecino(d));
			end
			%fprintf('6');
			for d=1:dim
				temp_1=T(coord_vec_vec(2*d-1,1),coord_vec_vec(2*d-1,2));
				temp_2=T(coord_vec_vec(2*d,1),coord_vec_vec(2*d,2));
				m(d)=min(temp_1,temp_2);
			end
			%fprintf('7');
			[m,ind_m]=sort(m);
			k=2;
			u(1)=m(1)+delta/F_curr;
			u(k)=u(1);
			me_pase_de_dim=0;
			% hay que revisar esta cuenta, da boniatos
			while u(k-1)>m(k) & ~me_pase_de_dim
				%fprintf('8');
				%u(k)=sum(m(1:k))+sqrt( sum(m(1:k))^2 + k*delta^2/F_curr + k*sum(m(1:k).^2) );
				u(k)=sum(m(1:k))+sqrt( sum(m(1:k))^2 + k*delta^2/F_curr - k*sum(m(1:k).^2) );
				u(k)=u(k)/k;
				if k<dim-2
					k=k+1;
				else
					me_pase_de_dim=1;
				end
			end
			%T(current_vecino(1),current_vecino(2))=u(k);
			% REVISAR: esto puede estar mal
			% estoy cambiando el valor aunque no sea NW tendria que revisar esto.
			% capaz que hay otro lugar mejor para hacerlo...
			% k o k-1???
			if tipos(current_vecino(1),current_vecino(2))==1
				T(current_vecino(1),current_vecino(2))=u(k);
			end
			%fprintf('9');
			% keyboard
		end
	else
		stop=1;
	end
	iter=iter+1;
	%fprintf('paso: %d \n',iter);
	if (plot)
		figure(1);
		mesh(T);
		axis([0 50 0 50 0 50])
		drawnow;
	end
	%colormap(gray(256));
	%     figure(2);
	%     mapa_clas = [1 1 1;hsv(3)];
	%     imshow(tipos+1,mapa_clas);
	%    h=gcf;
	%    set(h,'Position',[148 104 528 528]);
	%    hh=get(h,'Children')
	%    set(hh,'Position',[0.1 0.1 0.8 0.8]);
	%    pause;
end
if (plot)
	figure
	%image(255./(max(max(T))-min(min(T)))*(T-min(min(T))));
	%colormap(gray(256))
	mesh(T);
end