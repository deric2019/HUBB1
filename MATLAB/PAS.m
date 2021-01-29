function E2=PAS(E1,L,N,a,lambda_noll,n_medium)

% Varje sampelpunkt i k-planet motsvarar en plan våg med en viss riktning (kx,ky,kz)
delta_k= 2*pi/(N*a); % samplingsavstånd i k-planet *** Klar
kxvekt=-N/2*delta_k:delta_k:(N/2-1)*delta_k; % vektor med sampelpositioner i kx-led
kyvekt=kxvekt; % och ky-led
[kxmat,kymat]=meshgrid(kxvekt,kyvekt); % k-vektorns x- resp y-komponent i varje sampelpunkt i k-planet

k=2*pi./lambda_noll.*n_medium; % k-vektorns längd (skalär) för en plan våg i ett material med brytningsindex n_medium *** Ej klar
kzmat=sqrt(k.^2-kxmat.^2-kymat.^2); % k-vektorns z-komponent i varje sampelpunkt i k-planet *** Klar (Obs! Matlab tillåter att en skalär direkt adderas/subtraheras med matris, man behöver alltså tex inte skriva "skalär*ones(N,N)-matris")

fasfaktor_propagation=exp(1i.*kzmat.*L); % faktorn varje sampelpunkt i k-planet (som ju motsvarar plan våg i viss riktning) multas med för att propagera sträckan L i z-led *** Klar

A= a^2/(2*pi)^2 * fft2c(E1); % Planvågsspektrum i Plan 1 *** Klar

B=A.*fasfaktor_propagation; % Planvågsspektrum i Plan 2 (Planvågsspektrum i Plan 1 multat med fasfaktorn för propagation för varje plan våg)

E2= delta_k^2* N^2 .*ifft2c(B); % *** Klar
end
