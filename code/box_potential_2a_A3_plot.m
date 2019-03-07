function box_potential_2a_A3_plot(V0)

% clc;
% clear all;
% close all;

%Program Parameters
nx=200;
ny=100;
%V0=1;

%Variables
sigma_box=1e-2;
sigma_main=1;
cMap=ones(ny,nx);
B=zeros(1,nx*ny);
G=sparse(nx*ny,nx*ny);


%Set conductivity map
cMap(:,:)=sigma_main;
cMap(1:40,81:120)=sigma_box;
cMap(61:100,81:120)=sigma_box;

for i=1:nx
    for j=1:ny
        N=j+(i-1)*ny;
        
        if(i==1)%left
            
            %G(N,:)=0;
            G(N,N)=1;
            B(N)=V0;
            
        elseif(i==nx)%right
            
            %G(N,:)=0;
            G(N,N)=1;
            B(N)=0;
            
        elseif(j==1)%bottom
            Nyu=j+1+(i-1)*ny;
            Nxu=j+(i-1+1)*ny;
            Nxd=j+(i-1-1)*ny;
            
            Ryu=(cMap(j,i)+cMap(j+1,i))/2.0;
            Rxu=(cMap(j,i)+cMap(j,i+1))/2.0;
            Rxd=(cMap(j,i)+cMap(j,i-1))/2.0;
            
            %G(N,:)=0;
            G(N,N)=-(Ryu+Rxu+Rxd);
            G(N,Nyu)=Ryu;
            G(N,Nxu)=Rxu;
            G(N,Nxd)=Rxd;
            
            
        elseif(j==ny)%top
            Nyd=j-1+(i-1)*ny;
            Nxu=j+(i-1+1)*ny;
            Nxd=j+(i-1-1)*ny;
            
            Ryd=(cMap(j,i)+cMap(j-1,i))/2.0;
            Rxu=(cMap(j,i)+cMap(j,i+1))/2.0;
            Rxd=(cMap(j,i)+cMap(j,i-1))/2.0;
            
            %G(N,:)=0;
            G(N,N)=-(Ryd+Rxu+Rxd);
            G(N,Nyd)=Ryd;
            G(N,Nxu)=Rxu;
            G(N,Nxd)=Rxd;
            
        else%middle
            Nyu=j+1+(i-1)*ny;
            Nyd=j-1+(i-1)*ny;
            Nxu=j+(i-1+1)*ny;
            Nxd=j+(i-1-1)*ny;

            Ryu=(cMap(j,i)+cMap(j+1,i))/2.0;
            Ryd=(cMap(j,i)+cMap(j-1,i))/2.0;
            Rxu=(cMap(j,i)+cMap(j,i+1))/2.0;
            Rxd=(cMap(j,i)+cMap(j,i-1))/2.0;
            
            %G(N,:)=0;
            G(N,N)=-(Ryu+Ryd+Rxu+Rxd);
            G(N,Nyu)=Ryu;
            G(N,Nyd)=Ryd;
            G(N,Nxu)=Rxu;
            G(N,Nxd)=Rxd;
            
        end
            
            
            
        
    end
end



V=G\B';

for i=1:nx
    for j=1:ny


        N=j+(i-1)*ny;

        Vmat(j,i)=V(N);

    end
end

[Ex,Ey]=gradient(-Vmat);
Emag=sqrt(Ex.^2+Ey.^2);
Ex_scale=Ex.*1e9.*200./201;
Ey_scale=Ey.*1e9.*100./101;

% CurrentLeft=sum(cMap(:,1).*Ex(:,1));
% CurrentRight=sum(cMap(:,nx).*Ex(:,nx));
% 
% fprintf('The relative current througth the left contact is: %4.5f.\n',CurrentLeft);
% fprintf('The relative current througth the right contact is: %4.5f.\n',CurrentRight);

% startlabel=3;
% figure(startlabel+1);
% surf(cMap,'EdgeColor','flat');
% title('Figure 4: Conductivity Map');
% ylabel('y position');
% xlabel('x position');
% cb4=colorbar;
% cb4.Label.String = '\sigma (S/m)';
% view(0,90);


figure(5);
%surf(Vmat,'EdgeColor','flat');
contourf(Vmat,30);
title('Figure 5: Voltage Map with Bottle Neck');
ylabel('y position bin');
xlabel('x position bin');
cb6=colorbar;
cb6.Label.String = 'V (V)';
view(0,90);

figure(6);
quiver(Ex,Ey);
title('Figure 6: Electric Field Vector Map with Bottle Neck');
ylabel('y position bin');
xlabel('x position bin');
xlim([0 nx]);
ylim([0 ny]);

% figure(startlabel+4);
% surf(Emag);
% title('Electric Field Magnitude Map');
% ylabel('y position');
% xlabel('x position');

% figure(startlabel+4);
% quiver(cMap.*Ex,cMap.*Ey);
% title('Figure 7: Current Density Vector Map');
% ylabel('y position');
% xlabel('x position');
% xlim([0 nx]);
% ylim([0 ny]);

end


