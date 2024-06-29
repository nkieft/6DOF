function out = ENUToECEF(X_e, Y_e, Z_e, in)
%ENUTOECEF transforms a vector from the ENU frame to the ECI frame. Does
%not account for movement of the frames relative to one another. 

%Latitude and longitude
[lon, lat, h] = ECEFtoGeo(X_e, Y_e, Z_e);
lon = deg2rad(lon);
lat = deg2rad(lat);

%Transformation matrix
R = [-sin(lon), -cos(lon)*sin(lat), cos(lon)*cos(lat);
    cos(lon), -sin(lon)*sin(lat), sin(lon)*cos(lat)
    0, cos(lat), sin(lat)];

%Output vector
out = R*in;


end