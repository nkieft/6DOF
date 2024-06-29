function [lon, lat, h] = ECEFtoGeo(X, Y, Z)

a = 6378137.0;
b = 6356752.3142;
e2 = (a.^2-b.^2)./a.^2;
ep2 = (a.^2-b.^2)./b.^2;
p = sqrt(X.^2+Y.^2);
F = 54.*b.^2.*Z.^2;
G = p.^2+(1-e2).*Z.^2-e2.*(a.^2-b.^2);
c = e2.^2.*F.*p.^2./G.^3;
s = (1+c+sqrt(c.^2+2.*c)).^(1./3);
k = s+1+1./s;
P = F./3./k.^2./G.^2;
Q = sqrt(1+2.*e2.^2.*P);
r0 = -P.*e2.*p./(1+Q)+sqrt(1./2.*a.^2.*(1+1./Q)-P.*(1-e2).*Z.^2./Q./(1+Q)-1./2.*P.*p.^2);
U = sqrt((p-e2.*r0).^2+Z.^2);
V = sqrt((p-e2.*r0).^2+(1-e2).*Z.^2);
z0 = b.^2.*Z./a./V;
h = U.*(1-b.^2./a./V);
lat = atand((Z+ep2.*z0)./p);
lon = atan2d(Y, X);

end