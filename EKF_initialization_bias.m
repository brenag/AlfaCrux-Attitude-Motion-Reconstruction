% Script used to initialize the EKF covariance matrices based on selected
% parameters


%% Measurements error

    meas_sigma = 200; % sigma of the measurements error in nT
    b_sigma=meas_sigma/20000;  % relative error

%% Covariance Matrices and State Vector initialization

    sigma_d = 1e-10;                       %Diagonal of dynamics error matrix D
    sigma_dr = 1e-6;                      %Error in magnetic dipole
    sigma_db = 1e-4;                      %Error in magnetemeter bias
    sigma_q0 = 10^-4;                     %Diagonal of the initial covariance matrix P0
    sigma_w0 = 10^-4;                     %Diagonal of the initial covariance matrix P0
    sigma_r0 = 10^-4;                     %Diagonal of the initial covariance matrix P0
    sigma_mb0 = 10^-4;                     %Diagonal of the initial covariance matrix P0
    
    Kalman_parameters = cell(2,1);        %All the matrix are in the cell
    
    EKF_parameters{1,1} = diag([ones(3,1)*sigma_d^2; ones(3,1)*sigma_dr^2; ones(3,1)*sigma_db^2]);    %Matrix D
    EKF_parameters{2,1} = diag(ones(3,1)*b_sigma^2);                           %Matrix R
    P = diag([ones(3,1)*sigma_q0^2; ones(3,1)*sigma_w0^2; ones(3,1)*sigma_r0^2; ones(3,1)*sigma_mb0^2]);   %Matrix P0
    
 %% Proposed initial conditions while not using DE algorithm

    % State(1:4, 1) = [1; 0; 0; 0]; %Initial guess for the quaternion
    % State(5:7, 1) = [-1; -2; -0.5]*pi/180; %Initial guess for the angular velocity vector
    % State(8:10, 1) = [0; 0; 0]; %Initial guess for the magnetic dipole
