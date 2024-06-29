function dispersion = Post_disp(longF, latF, L_lat, L_lon, Apogee)
%Determines East/North position of each landing site relative to launch
%site. Looked it up at some point, needa find source for the math
for i = 1:length(latF)

    lat1 = deg2rad(L_lat);
    lat2 = deg2rad(latF(i));
    dlat = lat2-lat1;
    dlon = deg2rad(longF(i)-L_lon);

    a = sin(dlat/2)^2+cos(lat1)*cos(lat2)*sin(dlon/2)^2;
    c = 2*atan2(sqrt(a), sqrt(1-a));
    d(i) = 6378100*c;

    dN(i) = dlat*6378100;
    dE(i) = sign(dlon)*sqrt(d(i)^2-dN(i)^2);

end

%Mean and standard deviation landing position and 
MN = mean(dN);
ME = mean(dE);
dist = sqrt((dN-MN).^2+(dE-ME).^2);
MD = mean(dist);
stdDist = std(dist);
rad = MD+stdDist;

dispersion.dN = dN;
dispersion.dE = dE;
dispersion.MN = MN;
dispersion.ME = ME;
dispersion.rad = rad;
dispersion.meanA = mean(Apogee);
dispersion.stdA = std(Apogee);

end
%{
%Plot 1-3 sigma
theta = 0:0.01:2*pi;
x1 = rad*cos(theta)+ME;
y1 = rad*sin(theta)+MN;
x2 = 2*rad*cos(theta)+ME;
y2 = 2*rad*sin(theta)+MN;
x3 = 3*rad*cos(theta)+ME;
y3 = 3*rad*sin(theta)+MN;

%Nothing but more plotting
figure
hold on
plot(dE, dN, 'ro');
plot(x1, y1, x2, y2, x3, y3);
xlabel("East/West");
ylabel("North/South");
%}