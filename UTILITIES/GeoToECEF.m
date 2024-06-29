function [X, Y, Z] = GeoToECEF(lat, lon, h)

a = 6378137.0;
f = 1/298.257223563;
b = a*(1-f);

e2 = 1-b^2/a^2;
N = a/sqrt(1-e2*sind(lat)^2);

X = (N+h)*cosd(lat)*cosd(lon);
Y = (N+h)*cosd(lat)*sind(lon);
Z = ((1-f)^2*N+h)*sind(lat);

end