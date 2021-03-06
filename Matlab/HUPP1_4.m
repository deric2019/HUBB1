clear all
clc
clf
close all
%full_white_value=64; % �ldre matlabversion - detta v�rde plottas som vitt (max) med image-kommandot
full_white_value=255; % nyare matlabversion - prova denna om din plot verkar m�rk!

N=1024; % NxN �r matrisstorleken (rekommenderad storlek N=1024)
sidlaengd_Plan1=4e-3; % det samplade omr�dets storlek (i x- eller y-led) i Plan 1 (rekommenderad storlek 4 mm)
a=sidlaengd_Plan1/N; % samplingsavst�nd i Plan 1 (och Plan 2 eftersom vi anv�nder PAS)
L=1000e-3; % propagationsstr�cka (dvs avst�nd mellan Plan 1 och 2)

lambda_noll=633e-9; % vakuumv�gl�ngd f�r r�tt ljus fr�n en HeNe-laser
n_medium=1; % brytningsindex f�r medium mellan Plan 1 och 2
k=2*pi/lambda_noll*n_medium; % k-vektorns l�ngd *** Klar

xvekt=-N/2*a:a:(N/2-1)*a; % vektor med sampelpositioner i x-led
yvekt=xvekt; % och y-led
[xmat,ymat]=meshgrid(xvekt,yvekt); % koordinatmatriser med x- och y-v�rdet i varje sampelposition
rmat=sqrt(xmat.^2+ymat.^2); % avst�ndet till origo i varje sampelpunkt. Observera att alla operationer �r elementvisa!

%******* F�lt i Plan 1
f_lins=1000e-3; % fokall�ngd p� linsen f�re Plan 1
T_lins=exp(-1i*k*rmat.^2/(2*f_lins)); % Transmissionsfunktion f�r en lins (linsen �r TOK)

D_apertur=2e-3;
T_apertur=rmat<(D_apertur/2); % Transmissionsfunktion f�r en cirkul�r apertur ("pupill") 

omega_in=1e-3; % 1/e2-radie (f�r intensiteten, dvs 1/e-radie f�r amplituden) f�r infallande Gaussiskt f�lt
E_in_gauss=exp(-rmat.^2/omega_in^2); % Infallande f�lt: Gaussiskt med plana v�gfronter och normalinfall (dvs konstant fas, h�r=0)

E_in_konstant=ones(N,N); % Infallande f�lt: Plan v�g med normalt infall

E1_gauss=E_in_gauss.*T_lins; % F�ltet i Plan 1 (precis efter linsen) f�r gaussisk str�le 
E1_cirkular=E_in_konstant.* T_lins .* T_apertur; % F�ltet i Plan 1 (precis efter linsen) f�r konstant f�lt som passerat genom cirkul�r apertur *** Ej klar 
E1_in_hermitegauss = E_in_gauss.*xmat; 

% E1=E1_gauss; % V�lj fall!
E1=E1_cirkular; 
% E1=E1_in_hermitegauss; 


I1=abs(E1).^2; % intensiteten �r prop mot kvadraten p� f�ltets amplitud (normalt struntar man i proportionalitetskonstanten)
figure(1)
image(xvekt*1e3,yvekt*1e3,I1/max(max(I1))*full_white_value)
title('Intensitet i Plan 1. Verkar OK, eller?')
xlabel('x [mm]')
ylabel('y [mm]')
colormap(gray) 
drawnow
axis('equal')

figure(2)
imagesc(xvekt*1e3,yvekt*1e3,angle(E1))
title('Fas i Plan 1. Verkar OK, eller?')
xlabel('x [mm]')
ylabel('y [mm]')
colormap(gray) 
colorbar
drawnow
axis('equal')


%**** Och nu propagerar vi till Plan 2!

E2=PAS(E1,L,N,a,lambda_noll,n_medium); % Propagation med PAS-funktionen *** Klar

I2=abs(E2).^2; 
mattnadsfaktor_plot=90; % anger hur m�nga g�nger maxv�rdet ska vara m�ttat i plotten (>1, kan vara bra om man vill se svagare detaljer)

figure(3)
image(xvekt*1e3,yvekt*1e3,I2/max(max(I2))*full_white_value*mattnadsfaktor_plot)
title(['Intensitet efter ' num2str(L*1e3) ' mm propagation (mattnadsfaktor=' num2str(mattnadsfaktor_plot) '). Verkar OK, eller?'])
xlabel('x [mm]')
ylabel('y [mm]')
colormap(gray) 
drawnow
axis('equal')

figure(4)
plot(xvekt*1e3,I2(N/2+1,:))
title(['Intensitet l�ngs x-axeln efter ' num2str(L*1e3) ' mm propagation. Verkar OK, eller?'])
xlabel('x [mm]')
drawnow
omega = [0.03868e-3 0.392e-3];
D_spot = 2.*omega; % Enligt definition
C=D_spot.*D_apertur./(lambda_noll.*[100e-3 1000e-3]);


