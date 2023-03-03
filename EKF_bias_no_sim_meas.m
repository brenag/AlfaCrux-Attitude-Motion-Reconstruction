function [Est, P] = EKF_bias_no_sim_meas(Est, P, B_sat, B_eci, Kalman_parameters,h, index)
% Extended Kalman Filter for Attitude Motion Estimation 
% State Vector x = [q, omega, m_res] 


%% Kalman filter matrix initialization
D = Kalman_parameters{1,1};

meas_sigma = 200; 
b_sigma = meas_sigma/norm(B_eci);
R = diag(ones(3,1)*b_sigma^2); 
% R = Kalman_parameters{2,1};
% quat2dcm

global J worb
worb = 0.0011;

% Remove this variable if the bias wont affect the filter
real_time_bias = Est(11:13);
%disp(real_time_bias);


% Determine the IGRF model for the sattelite body reference frame
B_model = quatrotate(Est(1:4)', B_eci')';
B_body = B_model/ norm(B_model) + real_time_bias/norm(B_model);

%% Attitude motion equations integration
%Returns the State Vector and Quaternion representation of the attitude
%matrix 

[A, Est] = motion_integration_bias(Est, B_model, index);

%% Measurement and dynamic matrix calculation

W_B = [0     -B_model(3)  B_model(2);
     B_model(3)  0     -B_model(1);
    -B_model(2) B_model(1)  0];       % [B0 x]

W_b = [0     -B_body(3)  B_body(2);
     B_body(3)  0     -B_body(1);
    -B_body(2) B_body(1)  0];       % [B0 x]

W = [0     -Est(7)  Est(6);
   Est(7)  0     -Est(5);
  -Est(6) Est(5)  0];               % [Omega x]


% Matrix of linearized gyroscopic torque
Jw = J*Est(5:7);

W_Jw = [0    -Jw(3)  Jw(2);
      Jw(3)   0   -Jw(1);
     -Jw(2) Jw(1)    0];        % [Jw x]
 
F_gir = 2*(W_Jw*W - W*J*W);
% F_gir = -2*(W_Jw - W*J);

% Matrix of linearized gravity gradient torque
e3 = [0; 0; 1];

Ae3 = A*e3;
JAe3 = J*Ae3;

W_Ae3 = [0    -Ae3(3)  Ae3(2);
      Ae3(3)   0   -Ae3(1);
     -Ae3(2) Ae3(1)    0];        % [Ae3 x]

W_JAe3 = [0    -JAe3(3)  JAe3(2);
      JAe3(3)   0   -JAe3(1);
     -JAe3(2) JAe3(1)    0];        % [JAe3 x]

F_gr = 6*worb^2*(W_Ae3*J*W_Ae3 - W_JAe3*W_Ae3);


% Matrix of magnetic dipole torque  

W_m = [0     -Est(10)  Est(9);
       Est(10)  0     -Est(8);
      -Est(9) Est(8)  0];           %[m x]
  
% F=mxB=-Bxm
F_m = 2*W_m*W_B*1e-9;



%% Matrix of Dynamics
F = [-W     eye(3)*0.5     zeros(3,6);
    inv(J)*(F_gr + F_m)    inv(J)*F_gir    -inv(J)*W_B*1e-9    zeros(3,3);
    zeros(6,12)];
%-inv(J)*W_B*1e-9*0  

%% Matrix of Measurements (3x9)
H = [2*W_b zeros(3,3) zeros(3,3) eye(3)];


%% State Transition Matrix
%exp(F(t_k - t_{k-1}) Fi \approx I_n + delta_t * F
Fi = eye(12,12) + F*h;

%% Matrix of dynamical noise
G = [zeros(3,9);[inv(J) zeros(3,6)];[zeros(3,6) eye(3,3)];[zeros(3,6) eye(3,3)]];
% instead of the integral, multiplies by the delta_t (h)
Q = Fi*G*D*G'*Fi'*h;

%% Covariance matrix calculation
P = Fi*P*Fi'+Q;
K = P*H'/(H*P*H'+R);
P = P-K*H*P;

%% State vector correction

% Bias correction
% B_sat = B_sat - real_time_bias;

dX = K*(B_sat/norm(B_model)-B_body);

Est(5:10) = Est(5:10)+dX(4:9);
Est(11:13) = Est(11:13)+dX(10:12)*norm(B_model);
% Est(5:7) = Est(5:7)+dX(4:6);
Est(1:4) = quatnormalize(quatmultiply(Est(1:4)',[sqrt(1-sum(dX(1:3).^2)) dX(1:3)']))';

%% Observability test

%Ob = obsv(F,H);
%rank(Ob)
%disp(Ob);
%unobsv = length(F) - rank(Ob)

end
