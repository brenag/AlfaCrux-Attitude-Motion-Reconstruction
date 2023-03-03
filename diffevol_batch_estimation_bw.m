% Script used to initialize some Differential Evolution parameters and also
% generate the plots related to the cost function

global B_eci B_sat t J worb

% AlfaCrux Inertial Matrix
J = [0.00183508557, -0.00000456992, 0.00000877987;
    -0.00000456992, 0.00185284505, 0.00000153435;
    0.00000877987, 0.00000153435, 0.00184586679];

% Orbit local normal and orbit angular rate

worb = 0.0011; % 2*pi/T , T= 95min

%-------------------------------------------------
% Initial conditions
%q0 = [0.5; 0.5; 0.5; 0.5];
%w0 = [-1; -2; -0.5]*pi/180;


% State vector x = [quaternion; angular rate; dipole]
%x0 = [q0; w0; m0];

%-----------------------------------------------------------------------------
%% Differential Evolution Parameters 
%-----------------------------------------------------------------------------

F = 0.1; % Differential weight
CR = 0.9; % Crossover Probability
NP = 160; % Number of the population vectors
D = 16; % State dimension
G = 500; % number of generations

%-----------------------------------------------------------------------------
%% DE execution and plotting
%-----------------------------------------------------------------------------

x_est = differential_evolution_bw(@cost_func_bw,F, CR, NP, D, G);
cost_func_plot_bw(x_est)

