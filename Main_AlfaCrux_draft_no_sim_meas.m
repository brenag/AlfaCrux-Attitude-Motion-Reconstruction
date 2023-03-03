% Data proceccing of AlfaCrux satellite

clc
clear all
close all

global t J worb Omega_meas

%-----------------------------------------------------------------------------
%% Extract data from the .mat workspace
%-----------------------------------------------------------------------------


% Load the MATLAB workspace with the measurements
% Make sure to add the right path 
load('/home/brenag/Desktop/AlfaCrux-Bias-Test/test/data/August_7/Data_08_07_2022.mat')

% Remember to change the name of the file with the corrected measurements
load('/home/brenag/Desktop/AlfaCrux-Bias-Test/test/data/August_7/Mag_measurements_without_bias_08_07_2022.mat')



% Epoch - Time of the first data sample
year_meas = Meas(1,1);                    
month_meas = Meas(1,2);
day_meas = Meas(1,3);
Chass = Meas(1,4);
Mins = Meas(1,5);
Secs = Meas(1,6);

t_begin = [year_meas month_meas day_meas Chass Mins Secs];

% Convert the time of each sample to Julian date 
Time = jday(Meas(:,1), Meas(:,2), Meas(:,3), Meas(:,4), Meas(:,5), Meas(:,6)); 

% Calculate the time spent since epoch (first sample) from each sample (in seconds)
Time = (Time - Time(1)) * 24 * 60 * 60; 
t = Time;

% Get the magnetometer measurements (in mG) and convert to nT
B_sat = Meas(:,10:12) * 100; 

Omega_meas = [Meas(:,7) Meas(:,8) Meas(:,9)] * pi/180;
Omega_meas = Omega_meas';


%-----------------------------------------------------------------------------
%% Extract TLE from the .txt file
%-----------------------------------------------------------------------------

% Read the lines of the TLE file
lines = readlines("/home/brenag/Desktop/AlfaCrux-Bias-Test/test/data/August_7/TLE_08_07_2022.txt");
longstr1 = char(lines(1,1));
longstr2 = char(lines(2,1));
  

%-----------------------------------------------------------------------------
%% Initiate some variables
%-----------------------------------------------------------------------------

% Some flags for the SGP4, IGRF, and Sun direction models (embed_func)
flag1 = 1;                            
flag2 = 1;
flag3 = 1;  

% Set the length of the telemetry data array
N=length(t); 

S_orbital = zeros(3,N); % Sun direction in the orbital reference frame
B_orbital = zeros(3,N); % Earth magnetic field induction in the orbital reference frame
S_eci = zeros(3,N);     % Sun direction in the inertial reference frame
B_eci = zeros(3,N);     % Earth magnetic field induction in the inertial reference frame
xsat_eci = zeros(3,N);  % Radius-vector in ECI
vsat_eci = zeros(3,N);  % Velocity-vector in ECI
xsat_ecf = zeros(3,N);  % Radius-vector in ECF
vsat_ecf = zeros(3,N);  % Velocity-vector in ECF
Time=t;

% AlfaCrux Tensor of Inertia
J = [0.00183508557, -0.00000456992, 0.00000877987;
    -0.00000456992, 0.00185284505, 0.00000153435;
    0.00000877987, 0.00000153435, 0.00184586679];



% cycle for obtaining the magnetic field induction vector and Sun direction
% in ECI
for nom=1:1:length(t) 
    [omega_orbital, S_orbital(:,nom), B_orbital(:,nom), S_eci(:,nom), B_eci(:,nom), xsat_eci(:,nom), vsat_eci(:,nom), gst, SHADOW(:,nom), ~, ~, ~, xsat_ecf(:,nom), vsat_ecf(:,nom)] = embed_func(t(nom), t_begin, flag1, flag2, flag3, longstr1, longstr2);
end


% Orbital Velocity modulus 
worb = norm(omega_orbital);

%-----------------------------------------------------------------------------
%% Differential Evolution 
%-----------------------------------------------------------------------------

% The process of the differential evolution might take some time so the ';'
% is removed for the user to follow the progress of the execution in the
% command window

% Remove if you decide to use true magnetometer data instead of the one
% with the bias correction
%B_sat = B;

% The user can also use '%' to simulate with previous DE results

%diffevol_batch_estimation 
%x_est  

%-----------------------------------------------------------------------------
%% Kalman Filter
%-----------------------------------------------------------------------------


State = zeros(13,length(t)); 
%State(:,1) = x_est;

% Set the initial state with the previous differential evolution results
%State(:,1) = [0.211894384478143 %cost = 0.0448
%0.460403488710936
%0.781143145468144
%0.364615939993799
%-0.00100948247683568
%-0.00103850967645277
%-0.00525154252767399
%0.0113677209470236
%0.0235756628571977
%-0.00188608489393702
%21622.0110889414
%-5037.19279940888
%24712.9762964445];

State(:,1) = [0.388369471292978
0.0663635863422363
-0.900652072342432
0.183278129520932
0.00369259286700461
-0.000104476265342422
-0.00128799629168509
0.0188386585940047
0.00288699708811700
-0.0198021840529504
-3098.82417657773
12200.1940456180
437.155923421165];

State(1:10,1) = [0.2065
    0.9307
    0.2235
    0.2031
   -0.0104
   -0.0008
    0.0086
    0.0206
    0.0029
   -0.0222];

% Initial conditions estimated by constant bias
%State(11:13,1) = [18645.3631040990
%-12977.6404446235
%-28006.6829528688];



% Initialization of the EKF covariance matrices
EKF_initialization_bias;


% Extended Kalman filter iterative cycle 
for sample = 1:(length(t)-1)
    
    
    % Rotate the Magnetic Field vector by the current quaternion 
    %B_model(:,sample) = quatrotate(State(1:4, sample)', B_eci(:,sample)');

    % Once we've got telemetry with high sample rate, integrate the motion
    % to compare the results with the filter
    
    % Estimate the current state using the EKF
    [State(:, sample + 1), P] = EKF_bias_no_sim_meas(State(:, sample), P, B_sat(sample, :)', B_eci(:,sample), EKF_parameters, 1, sample);

end 


%-----------------------------------------------------------------------------
%% Telemetry and Orbit Plots
%-----------------------------------------------------------------------------


% Earth orbit plot 

figure('Color',[1 1 1])
plot3(xsat_eci(1,:),xsat_eci(2,:),xsat_eci(3,:))
axis equal
hold on
plot3(xsat_eci(1,1),xsat_eci(2,1),xsat_eci(3,1),'or')
grid on

load('topo.mat','topo','topomap1'); % For Earth plotting
whos topo topomap1;
[x,y,z] = sphere(50);
R_earth = 6.371e6; %
% [x, y, z] = EarthRotation(x, y, z, i*dt);
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;
surface(x*R_earth,y*R_earth,z*R_earth,props);
light('position',[-1 0 1]);
light('position',[-1.5 0.5 -0.5], 'color', [.6 .2 .2]);
axis equal
view(45,45)
xlabel('X, m')
ylabel('Y, m')
zlabel('Z, m')


% Altitude

figure('Color',[1 1 1])
plot((xsat_eci(1,:).^2+xsat_eci(2,:).^2+xsat_eci(3,:).^2).^(1/2)-6.4e6,'LineWidth',2)
xlabel('Time, s')
ylabel('Altitude, m')
grid on


% Magnetic Field in ECI frame

figure('Color',[1 1 1])
plot(Time,B_eci(1,:),'b','LineWidth',2);
hold on
plot(Time,B_eci(2,:),'g','LineWidth',2);
plot(Time,B_eci(3,:),'r','LineWidth',2);
xlabel('Time, s')
ylabel('Magnetic field in ECI, nT')
grid on


% Magnetometer measurements with bias correction

%load('Mag_measurements without bias.mat')

figure('Color',[1 1 1])
plot(Time,B(:,1),'b','LineWidth',2);
hold on
plot(Time,B(:,2),'g','LineWidth',2);
plot(Time,B(:,3),'r','LineWidth',2);
xlabel('Time, s')
ylabel('Magnetometer measurements without bias, nT')
grid on
legend('B_x','B_y','B_z')


% Magnetometer measurements

figure('Color',[1 1 1])
plot(Time,B_sat(:,1),'b','LineWidth',2);
hold on
plot(Time,B_sat(:,2),'g','LineWidth',2);
plot(Time,B_sat(:,3),'r','LineWidth',2);
xlabel('Time, s')
ylabel('Magnetometer measurements, nT')
grid on
legend('B_x','B_y','B_z')


% Comparison between the modulus of the magnetic field in ECI and the
% measurements with and without bias

figure('Color',[1 1 1])
plot(Time,(B_eci(1,:).^2+B_eci(2,:).^2+B_eci(3,:).^2).^(0.5),'b','LineWidth',2);
hold on
plot(Time,(B_sat(:,1).^2+B_sat(:,2).^2+B_sat(:,3).^2).^(0.5),'*-r','LineWidth',2);
plot(Time,(B(:,1).^2+B(:,2).^2+B(:,3).^2).^(0.5),'*-g','LineWidth',2);
xlabel('Time, s')
ylabel('Magnetic field value, nT')
legend('Model','Measurements','Measurements without bias')
grid on


% Sun direction vector

figure('Color',[1 1 1])
plot(Time,S_eci(1,:),'b','LineWidth',2);
hold on
plot(Time,S_eci(2,:),'g','LineWidth',2);
plot(Time,S_eci(3,:),'r','LineWidth',2);
xlabel('Time, s')
ylabel('Sun direction vector')
grid on


% Gyroscope measurements (in deg/s)

Meas(:,7) = detrend(Meas(:,7));
Meas(:,8) = detrend(Meas(:,8));
Meas(:,9) = detrend(Meas(:,9));

figure('Color',[1 1 1])
plot(Time,Meas(:,7),'b','LineWidth',2);
hold on
plot(Time,Meas(:,8),'g','LineWidth',2);
plot(Time,Meas(:,9),'r','LineWidth',2);
xlabel('Time, s')
ylabel('Angular velocity measurements, deg/s')
grid on


% Angular momentum (in kg*m²/s)

for k = 1:1:length(Meas(:,7))
    K(:,k) = J * Meas(k,7:9)' * pi/180;
    K_norm(k) = norm(K(:,k));
end

figure('Color',[1 1 1])
plot(Time,K(1,:),'b','LineWidth',2);
hold on
plot(Time,K(2,:),'g','LineWidth',2);
plot(Time,K(3,:),'r','LineWidth',2);
plot(Time,K_norm,'k','LineWidth',2);
xlabel('Time, s')
ylabel('Angular momentum, kg*m^2/s')
grid on


% Current from each solar panel

figure('Color',[1 1 1])
plot(Time,Meas(:,13),'b','LineWidth',2);
hold on
plot(Time,Meas(:,14),'g','LineWidth',2);
plot(Time,Meas(:,15),'r','LineWidth',2);
xlabel('Time, s')
ylabel('Solar panels curents, mA')
grid on
legend('X','Y','Z')


% Voltage from each solar panel

figure('Color',[1 1 1])
plot(Time,Meas(:,16),'b','LineWidth',2);
hold on
plot(Time,Meas(:,17),'g','LineWidth',2);
plot(Time,Meas(:,18),'r','LineWidth',2);
xlabel('Time, s')
ylabel('Solar panels voltage, mV')
legend('X','Y','Z')
grid on


% Temperature of each solar panel

figure('Color',[1 1 1])
plot(Time,Meas(:,19),'b','LineWidth',2);
hold on
plot(Time,Meas(:,20),'--b','LineWidth',2);
plot(Time,Meas(:,21),'g','LineWidth',2);
plot(Time,Meas(:,22),'--g','LineWidth',2);
plot(Time,Meas(:,23),'r','LineWidth',2);
plot(Time,Meas(:,24),'--r','LineWidth',2);
xlabel('Time, s')
ylabel('Temperature, C')
legend('X+','X-','Y+','Y-','Z+','Z-')
grid on

%-----------------------------------------------------------------------------
%% Reconstructed Attitude Motion Plots
%-----------------------------------------------------------------------------

% EKF Magnetometer Bias

figure('Color',[1 1 1])
plot(t,State(11,:),'b','LineWidth',2);
hold on
plot(t,State(12,:),'g','LineWidth',2);
plot(t,State(13,:),'r','LineWidth',2);
xlabel('Time, s')
ylabel('Estimated magnetometer bias, nT')
grid on


% EKF Dipole Moment

figure('Color',[1 1 1])
plot(t,State(8,:),'b','LineWidth',2);
hold on
plot(t,State(9,:),'g','LineWidth',2);
plot(t,State(10,:),'r','LineWidth',2);
xlabel('Time, s')
ylabel('Estimated residual magnetic dipole, Am^2')
grid on


% EKF Angular Velocity

figure('Color',[1 1 1])
plot(t,State(5,:)*180/pi,'b','LineWidth',2);
hold on
plot(t,State(6,:)*180/pi,'g','LineWidth',2);
plot(t,State(7,:)*180/pi,'r','LineWidth',2);
xlabel('Time, s')
ylabel('Estimated Angular Velocity, deg/s')
grid on
legend('\omega_x','\omega_y','\omega_z')


% EKF Quaternion

figure('Color',[1 1 1])
plot(t,State(1,:),'b','LineWidth',2);
hold on
plot(t,State(2,:),'g','LineWidth',2);
plot(t,State(3,:),'r','LineWidth',2);
plot(t,State(4,:),'y','LineWidth',2);
xlabel('Time, s')
ylabel('Estimated Quaternion')
grid on
legend('q_0','q_1','q_2','q_3')



%% Bias and scale factor estimation

amp_DE = zeros(3,1);
amp_DE(1) = mean(State(5,:));
amp_DE(2) = mean(State(6,:));
amp_DE(3) = mean(State(7,:));
disp(amp_DE);

amp_meas = zeros(3,1);
amp_meas(1) = mean(Meas(:,7));
amp_meas(2) = mean(Meas(:,8));
amp_meas(3) = mean(Meas(:,9));
disp(amp_meas);

bias = zeros(3,1);
bias = amp_meas - amp_DE;
disp(bias);



scale_factor = zeros(3,1);

scale_factor(1) = (max(abs(State(5,:))) - abs(amp_DE(1))) / ( max(abs(Meas(:,7))) - abs(amp_meas(1)));
scale_factor(2) = (max(abs(State(6,:))) - abs(amp_DE(2))) / ( max(abs(Meas(:,8))) - abs(amp_meas(2)));
scale_factor(3) = (max(abs(State(7,:))) - abs(amp_DE(3))) / ( max(abs(Meas(:,9))) - abs(amp_meas(3)));

disp(scale_factor);

ang_vel = zeros(3, length(t));
ang_vel(1,:) = (Meas(:,7)' - bias(1))*scale_factor(1);
ang_vel(2,:) = (Meas(:,8)' - bias(2))*scale_factor(2);
ang_vel(3,:) = (Meas(:,9)' - bias(3))*scale_factor(3);


% Angular velocity without bias

figure('Color',[1 1 1])
plot(t,ang_vel(1,:)*180/pi,'b','LineWidth',2);
hold on
plot(t,ang_vel(2,:)*180/pi,'g','LineWidth',2);
plot(t,ang_vel(3,:)*180/pi,'r','LineWidth',2);
xlabel('Time, s')
ylabel('Angular Velocity without bias and scale factor, deg/s')
grid on
legend('\omega_x','\omega_y','\omega_z')

