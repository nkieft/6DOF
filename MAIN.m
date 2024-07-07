close all; clear; clc;

%% 6 Degree of Freedom Flight Simulator
%--------------------------------------------------------------------------
% Author: Nicholas Kieft
%
% Key References:
% Topics in Advanced Model Rocketry
% Openrocket Technical Documentation
% RASAero II Documentation
% Ender Kerr 6DOF Documentation
%--------------------------------------------------------------------------

%% Launch Site Inputs
%--------------------------------------------------------------------------
Launchsite_Longitude = -117.808244; %Longitude of launchsite [deg]
Launchsite_Latitude = 35.347142; %Latitude of launchsite [deg]
Launchsite_Elevation = 2080 / 3.281; %Elevation of launchsite [ft]

%% Launcher Inputs
%--------------------------------------------------------------------------
Launcher_Azimuth = 90; %Azimuth of launcher pivot axis [deg]
Launcher_Angle = 0; %Angle of launcher with respect to vertical [deg]
Rail_Length = 20; %Length of launch rail / tower [ft]
Rail_Contact_Force = 50; %Rocket-Rail Contact force
Rail_Friction_Coeff = 0; %Rocket-Rail Interface friction coefficient
%--------------------------------------------------------------------------

%% Environmental Inputs
%--------------------------------------------------------------------------
Manual_Winds = false; %If true, uses inputs below. If false, uses matlab HWM
%Wind Alt = [0 3 6 9 12 15 18 21 24 27 30 36 42 48] FOR REFERENCE
Winds_Mag = [3 3 6 9 8 5 3 2 4 3 12 24 37 45]; %Speed in knots
Winds_Dir = [100 100 130 140 140 150 240 230 150 270 290 290 290 290];%360 minus azimuth
Launch_Day = 110; %Day of the year of launch, necessary for HWM only
%--------------------------------------------------------------------------

%% Motor Inputs
%--------------------------------------------------------------------------
%Motor Inputs
Motor_File_Name = "MAXMAX.txt"; %Name of motor file to be used
Thrust_Angle = 0;
Thrust_Orientation = 0;
%--------------------------------------------------------------------------

%% Vehicle Inputs
%--------------------------------------------------------------------------
%Body
Nosecone_Length = 32.5; %Length of vehicle nosecone
Length_Vehicle = 144; %Length of entire Vehicle [in]
Radius_Vehicle = 3.25; %Radius of nosecone base [in]

%Fins
Fin_Root_Chord = 17; %Fin root chord length [in]
Fin_Tip_Chord = 6; %Fin Tip chord length [in]
Fin_Span = 5; %Fin semi-span [in]
Fin_Sweep = 11;
Fin_Cant = 0.0; %Fin cant angle [deg]
Fin_Thickness = 0.25; %Fin Thickness

%Coupler
Coupler_Location = 35; %Coupler distance from tip [in]
Coupler_Length = 15; %Length of empty volume containing parachute [in]
%--------------------------------------------------------------------------

%% Mass Inputs
%Mass and Inertia Overrides (Set to zero for no ovveride)
Initial_Mass_Override = 67.207*2.207;
Final_Mass_Override = 31.399*2.207;
Initial_Inertia_Override = 45.95;
Final_Inertia_Override = 33.57;
Initial_CG_Override = 88.18;
Final_CG_Override = 78.56;

%Normal Inputs
%--------------------------------------------------------------------------
% Initial Mass Distribution Inputs (Add all components contributing to mass)
% Format: s{i, :} = [Start X, End X, Start Diameter, End Diameter, Mass]
m_data_initial = table(); %Don't Remove this line
m_data_initial{1, :} = [0, 8, 0, 1.725, 2.51]; %Tip
m_data_initial{2, :} = [7.75, 32.25, 1.725, 6.5, 4]; %Nosecone
m_data_initial{3, :} = [32.5, 44.5, 6.5, 6.5, 4]; %Shoulder
m_data_initial{4, :} = [44.5, 48.5, 6.5, 6.5, 15.6/16]; %Switchband
m_data_initial{5, :} = [38, 55.125, 6.5, 6.5, 10.12]; %AV Bay
m_data_initial{6, :} = [48.5, 66.5, 6.5, 6.5, 4.48]; %Recovery Airframe
m_data_initial{7, :} = [63, 70.375, 6.5, 6.5, 6.8]; %Motor Bulkhead
m_data_initial{8, :} = [66.5, 137, 6.5, 6.5, 18.39]; %Motor Casing
m_data_initial{9, :} = [70, 81.5, 6, 6, 14.7]; %Grain 1
m_data_initial{10, :} = [81.5, 93, 6, 6, 14.7]; %Grain 2
m_data_initial{11, :} = [93, 104.5, 6, 6, 14.7]; %Grain 3
m_data_initial{12, :} = [104.5, 116, 6, 6, 14]; %Grain 4
m_data_initial{13, :} = [116, 127.5, 6, 6, 13.1]; %Grain 5
m_data_initial{14, :} = [127.5, 130, 6, 6, 11.8]; %Grain 6
m_data_initial{15, :} = [130, 137, 6, 6, 5.72]; %Nozzle
m_data_initial{16, :} = [144-17, 137, 6.5, 11, 7.16]; %Fins

% Burnout Mass Distribution Inputs (Input Any Components Contributing to Mass)
% Format: s{i, :} = [Start X, End X, Start Diameter, End Diameter, Mass]
m_data_final = table(); %Don't Remove this line
m_data_final{1, :} = [0, 8, 0, 1.725, 2.51]; %Tip
m_data_final{2, :} = [7.75, 32.25, 1.725, 6.5, 4]; %Nosecone
m_data_final{3, :} = [32.5, 44.5, 6.5, 6.5, 4]; %Shoulder
m_data_final{4, :} = [44.5, 48.5, 6.5, 6.5, 15.6/16]; %Switchband
m_data_final{5, :} = [38, 55.125, 6.5, 6.5, 10.12]; %AV Bay
m_data_final{6, :} = [48.5, 66.5, 6.5, 6.5, 4.48]; %Recovery Airframe
m_data_final{7, :} = [63, 70.375, 6.5, 6.5, 6.8]; %Motor Bulkhead
m_data_final{8, :} = [66.5, 137, 6.5, 6.5, 18.39]; %Motor Casing
m_data_final{9, :} = [70, 81.5, 6, 6, 0.5]; %Grain 1
m_data_final{10, :} = [81.5, 93, 6, 6, 0.5]; %Grain 2
m_data_final{11, :} = [93, 104.5, 6, 6, 0.5]; %Grain 3
m_data_final{12, :} = [104.5, 116, 6, 6, 0.5]; %Grain 4
m_data_final{13, :} = [116, 127.5, 6, 6, 0.5]; %Grain 5
m_data_final{14, :} = [127.5, 137, 6, 6, 0.5]; %Grain 6
m_data_final{15, :} = [130, 137, 6, 6, 5.72]; %Nozzle
m_data_final{16, :} = [144-17, 137, 6.5, 11, 7.16]; %Fins
%--------------------------------------------------------------------------

%% Recovery Inputs
%--------------------------------------------------------------------------
Ballistic = false; %TRUE for ballistic, FALSE for chuted recovery
Diameter_Drogue = 30; %Drogue Parachute Diameter [in]
Cd_Drogue = 2.2; %Drogue drag coefficient
Diameter_Main = 120; %Main Parachute Diameter [in]
Cd_Main = 2.2; %Main drag coefficient
Main_Deploy_Height = 1000; %Main deploy height [m]
%--------------------------------------------------------------------------

%% Simulation Inputs
%--------------------------------------------------------------------------
Termination_Time = 800; %Time at which simulation is forced to stop [s]
Run_Dispersion = false; %FALSE for single simulation, TRUE for monte carlo
Iteration_QTY = 1; %1 for single simulation, iteration qty for monte carlo
outfile = "Mamba"; %Name of output file
Max_Mach = 8; %Use 1.5x expected flight peak mach. Can do more, but will increase computation time
Max_Apogee = 150000 / 3.281; %Use 1.5x expected flight apogee. Can do more, but will increase computation time
%--------------------------------------------------------------------------

%% Monte Carlo Inputs
%--------------------------------------------------------------------------
Thrust_Sigma = 5; %Thrust std [%]
Mass_Sigma = 3; %Initial final mass std [%]
Launcher_Angle_Sigma = 1.5; %Launcher angle std  [deg]
Launcher_Azimuth_Sigma = 0; %Launcher angle std [deg]
Wind_Scale_Sigma = 20; %Wind scale std [%]
Wind_Offset_Sigma = 3; %Wind offset std [m/s]
C_Sigma = 10; %Aero coefficient std [%]
Cant_Sigma = 0.1; %Fin cant std [deg]
Thrust_Angle_Sigma = 0.1; %Off-axis angle [deg]
Thrust_Orientation_Sigma_Range = 360; %Off-axis orientation [deg]
%--------------------------------------------------------------------------

%% END OF INPUTS

%% Preliminary Calculations 
%--------------------------------------------------------------------------
%Add all subfolders to matlab path
addpath(genpath(pwd));

%Format Mass Data
m_data_initial = renamevars(m_data_initial, ["Var1", "Var2", "Var3", "Var4", "Var5"], ["X0", "XF", "D0", "DF", "M"]);
m_data_final = renamevars(m_data_final, ["Var1", "Var2", "Var3", "Var4", "Var5"], ["X0", "XF", "D0", "DF", "M"]);
m_data_initial{:, 1:4} = m_data_initial{:, 1:4}/39.37;
m_data_initial{:, 5} = m_data_initial{:, 5}/2.207;
m_data_final{:, 1:4} = m_data_final{:, 1:4}/39.37;
m_data_final{:, 5} = m_data_final{:, 5}/2.207;
Initial_Mass_Override = Initial_Mass_Override/2.207;
Final_Mass_Override = Final_Mass_Override/2.207;
Initial_CG_Override = Initial_CG_Override/39.37;
Final_CG_Override = Final_CG_Override/39.37;

%Calculate Recovery Areas
Area_Drogue = (Diameter_Drogue/2/39.37)^2*pi;
Area_Main = (Diameter_Main/2/39.37)^2*pi;

%Mass / Inertia Calculations
[CG_Initial, I_Initial, Mass_Initial, ~] = MassDistribution(m_data_initial, Length_Vehicle/39.37/1000, 1000, linspace(0, Length_Vehicle/39.37, 1000), Initial_Mass_Override, Initial_Inertia_Override, Initial_CG_Override);
[CG_Final, I_Final, Mass_Final, dmdx_f] = MassDistribution(m_data_final, Length_Vehicle/39.37/1000, 1000, linspace(0, Length_Vehicle/39.37, 1000), Final_Mass_Override, Final_Inertia_Override, Final_CG_Override);

%Load in lookup tables
F_T = readmatrix(Motor_File_Name); %Read thrust into lookup table
[C_b, C_c, Cp] = AeroImport(Max_Mach); %Import aero data from RASAero II
W_h = WindModel(Launchsite_Latitude, Launchsite_Longitude, Launch_Day, Winds_Mag, Winds_Dir, Manual_Winds, Max_Apogee); %Import wind model data

%Preload atmospheric data
h = 1:50:Max_Apogee;
[~, a, ~, rho_a] = atmosisa(h);
atm_data = [h', a', rho_a'];

%Aerodynamic Chord Calcs
ymac = Fin_Span/3*(Fin_Root_Chord+2*Fin_Tip_Chord)/(Fin_Tip_Chord+Fin_Root_Chord); %Calculate aero chord
X_f = Fin_Sweep/3*(Fin_Root_Chord+2*Fin_Tip_Chord)/(Fin_Root_Chord+Fin_Tip_Chord)+1/6*(Fin_Root_Chord^2+Fin_Tip_Chord^2+Fin_Root_Chord*Fin_Tip_Chord)/(Fin_Root_Chord+Fin_Tip_Chord);

%Split Thrust into time and force
F_T_t = F_T(:, 1);
F_T_T = F_T(:, 2);

%Vectorize standard deviations
if(Run_Dispersion)
    sigma = [Thrust_Sigma, Mass_Sigma, Mass_Sigma, Launcher_Azimuth_Sigma, Launcher_Angle_Sigma, Wind_Scale_Sigma, Wind_Offset_Sigma, C_Sigma, C_Sigma, Cant_Sigma, Thrust_Angle_Sigma, Thrust_Orientation_Sigma_Range];
else
    sigma = zeros(1, 12);
end

%% Dispersion Section
%--------------------------------------------------------------------------
for i = 1:Iteration_QTY %Loops for each iteration

    %List Simulation Number
    disp(strcat("SIMULATION # : ", string(i)));
    tic

    %Vary Monte Carlo Parameters
    [F_T_e, Mass_Initial_e, Mass_Final_e, Launcher_Azimuth_e, Launcher_Angle_e, C_b_e, C_c_e, Fin_Cant_e, Thrust_a_e, Thrust_o_e] = MonteCarlo(F_T_T, Mass_Initial, Mass_Final, Launcher_Azimuth, Launcher_Angle, C_b, C_c, Fin_Cant, Thrust_Angle, Thrust_Orientation, sigma);
    
    %Randomize wind
    W_r = WindUncertainty(W_h, Wind_Scale_Sigma, Wind_Offset_Sigma);

    %Calculate total motor impulse
    F_T_e(:, 2) = F_T_e;
    F_T_e(:, 1) = F_T_t;
    I_T = trapz(F_T_e(:, 2))*0.03;

    %Below 3 are constants to be passed into ode45 full of required data
    Vehicle_Data = FillVehicleData(Radius_Vehicle, Length_Vehicle, CG_Initial, CG_Final, Mass_Initial_e, Mass_Final_e, I_T, Fin_Root_Chord, Nosecone_Length, ymac, Fin_Cant_e, Main_Deploy_Height, Area_Main, Cd_Main, Area_Drogue, Cd_Drogue, I_Initial, I_Final, Thrust_a_e, Thrust_o_e,Fin_Tip_Chord, Fin_Sweep,Fin_Thickness, Fin_Span, X_f);
    Launchsite_Data = FillLaunchSiteData(Launchsite_Latitude, Launchsite_Longitude, Launcher_Azimuth_e, Launcher_Angle_e, Rail_Length, Launchsite_Elevation, Rail_Contact_Force, Rail_Friction_Coeff);
    Physical_Const = FillConst();

    %Initial conditions for ascent & termination conditions
    init = FillInit(Vehicle_Data, Launchsite_Data);
    OptA = odeset("Events", @(t, s)ApogeeEvent(t, s, Launchsite_Data, Ballistic), "RelTol", 1E-6, "AbsTol", 1E-6);

    %Ascent ODE loop
    [t, s] = ode45(@(t, s)Ascent(t, s, Vehicle_Data, Launchsite_Data, Physical_Const, F_T_e, C_b_e, C_c_e, Cp, W_r, atm_data), 0:0.01:Termination_Time, init, OptA);

    %Initial conditions for descent % termination conditions
    OptD = odeset("Events", @(t, s)LandingEvent(t, s, Launchsite_Data));
    init2 = [s(end, 11), s(end, 12), s(end, 13), s(end, 14), s(end, 15), s(end, 16), s(end, 18), s(end, 19), s(end, 20), s(end, 21), s(end, 22), s(end, 23)];

    %Chuted recovery only runs if ballistic recovery is not selected.
    if(~Ballistic)
        [t_p, s_p] = ode45(@(t, s)Descent(t, s, Vehicle_Data, Launchsite_Data, Physical_Const, W_r), [0 10000], init2, OptD);
    else
        t_p = t(end);
        s_p = init2;
    end

    %Save position data for dispersion plotting
    XD = [s(:, 18); s_p(:, 7)];
    YD = [s(:, 19); s_p(:, 8)];
    ZD = [s(:, 20); s_p(:, 9)];

    %Final longs/lats
    [lonD(i), latD(i), ~] = ECEFtoGeo(XD(end), YD(end), ZD(end));
    [~, ~, hD] = ECEFtoGeo(XD, YD, ZD);

    %Apogee for dispersion analysis
    Apogee(i) = max(hD);
    toc
end
%--------------------------------------------------------------------------

%% Mean Simulation 
%--------------------------------------------------------------------------
%Mean Impulse
I_T = trapz(F_T(:, 2))*0.03;

%Below 3 are constants to be passed into ode45 full of required data
Vehicle_Data = FillVehicleData(Radius_Vehicle, Length_Vehicle, CG_Initial, CG_Final, Mass_Initial, Mass_Final, I_T, Fin_Root_Chord, Nosecone_Length, ymac, Fin_Cant, Main_Deploy_Height, Area_Main, Cd_Main, Area_Drogue, Cd_Drogue, I_Initial, I_Final, Thrust_Angle, Thrust_Orientation, Fin_Tip_Chord, Fin_Sweep,Fin_Thickness, Fin_Span, X_f);
Launchsite_Data = FillLaunchSiteData(Launchsite_Latitude, Launchsite_Longitude, Launcher_Azimuth, Launcher_Angle, Rail_Length, Launchsite_Elevation, Rail_Contact_Force, Rail_Friction_Coeff);
Physical_Const = FillConst();

%Initial conditions for ascent & termination conditions
init = FillInit(Vehicle_Data, Launchsite_Data);
OptA = odeset("Events", @(t, s)ApogeeEvent(t, s, Launchsite_Data, Ballistic));

%Ascent ODE loop
[t, s] = ode45(@(t, s)Ascent(t, s, Vehicle_Data, Launchsite_Data, Physical_Const, F_T, C_b, C_c, Cp, W_h, atm_data), 0:0.01:Termination_Time, init, OptA);

%Initial conditions for descent % termination conditions
OptD = odeset("Events", @(t, s)LandingEvent(t, s, Launchsite_Data));
init2 = [s(end, 11), s(end, 12), s(end, 13), s(end, 14), s(end, 15), s(end, 16), s(end, 18), s(end, 19), s(end, 20), s(end, 21), s(end, 22), s(end, 23)];

%Chuted recovery only runs if ballistic recovery is not selected.
if(~Ballistic)
    [t_p, s_p] = ode45(@(t, s)Descent(t, s, Vehicle_Data, Launchsite_Data, Physical_Const, W_h), [0 10000], init2, OptD);
else
    %Set parachute outputs to final ascent for ballistic, makes
    %postprocessing easier
    t_p = t(end);
    s_p = init2;
end

%Plotting Extraction
XM = [s(:, 18); s_p(:, 7)];
YM = [s(:, 19); s_p(:, 8)];
ZM = [s(:, 20); s_p(:, 9)];

%Final longs/lats
[lonM, latM, hM] = ECEFtoGeo(XM, YM, ZM);

%Trajectory Plotting on Globe
%if = uifigure;
%earth = geoglobe(uif);
%hold(earth, 'on')
%geoplot3(earth, latM, lonM, hM, 'r', 'LineWidth', 1)

%% Sim Intermediate Value Extraction
%--------------------------------------------------------------------------
%Ascent analysis
for i = 1:height(s)
    output_a(i, :) = Post_A(t(i), s(i, :), Vehicle_Data, Launchsite_Data, Physical_Const, F_T, C_b, C_c, Cp, W_h, atm_data);
end

%Descent analysis
for i = 1:height(s_p)
    output_d(i, :) = Post_D(t_p(i, :), s_p(i, :), Vehicle_Data, Launchsite_Data, Physical_Const, W_h);
end

%Dispersion analysis
dispersion = Post_disp(lonD, latD, Launchsite_Latitude, Launchsite_Longitude, Apogee);
%--------------------------------------------------------------------------

%% App Prep and Formatting
%--------------------------------------------------------------------------
%Recovery Calculations
[Dsep.V_D, Dsep.V_M, Dsep.QTYPins, Dsep.W_BP, Dsep.F_sep] = RecoveryCalcs(Vehicle_Data, Launchsite_Data, dmdx_f, Coupler_Location, max(output_a(:, 9)), max(output_a(:, 20)), Coupler_Length, Mass_Final/2.207);

%Structural Calculations
[Struct.x, Struct.N, Struct.V, Struct.M] = StructuralCalcs(Vehicle_Data, max(output_a(:, 9)), max(output_a(:, 20)), max(F_T(:, 2)), Mass_Final/2.207, I_Final, CG_Final/39.37, dmdx_f);

%Dynamics Calculations
[Dyn.C1, Dyn.C2, Dyn.omega_n, Dyn.zeta, Dyn.omega_n_c, Dyn.zeta_c, Dyn.AR_res, Dyn.P, Dyn.t] = DynamicsCalcs(hM(1:length(s)), output_a(:, 18), Vehicle_Data, output_a(:, 24), output_a(:, 25), output_a(:, 26), output_a(:, 27), output_a(:, 28), output_a(:, 29), output_a(:, 3), output_a(:, 4), t, s(:, 4));

%Formatting table for custom app
output_d(:, 1) = output_d(:, 1)+t(end);
output = [output_a;output_d];
output = array2table(output);
output = renamevars(output, ["output1", "output2", "output3", "output4", "output5", "output6" ...
    "output7", "output8", "output9", "output10", "output11", "output12", "output13", "output14" ...
    "output15", "output16", "output17", "output18", "output19", "output20", "output21", "output22" ...
    "output23"], ...
    ["Time", "AoA", "CG", "CP", "Stability Margin", "Normal Force", "Drag Force", "Thrust Force", ...
    "Mach", "Latitude", "Longitude", "Mass", "Altitude", "q0", "q1", "q2", "q3", "Speed", "Acceleration" ...
    "Dynamic Pressure", "Roll Rate", "Pitch Rate", "Yaw Rate"]);

save(outfile, "output", "dispersion", "Dsep", "Struct", "Dyn");

%Removes Path (So you can modify folders etc)
rmpath(genpath(pwd));

%Boots up app
run("POSTPROCESS.mlapp");
