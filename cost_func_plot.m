function [F] = cost_func_plot(qw)


global B_eci B_sat xsat_eci t J worb Omega_meas  

e3 = xsat_eci(:,1)/norm(xsat_eci(:,1)); % Local Vertical ECI Vector

% Notice that the 'for loop' responsible for the integration of the attitude
% motion starts in index 2, so we must determine the magnetic field
% reference for the initial time


%Determine the IGRF model for the sattelite initial quaternion guess
B_model(:, 1) = quatrotate(qw(1:4)', B_eci(:,1)')' + qw(11:13);
B_model_norm(:, 1) = B_model(:, 1)/ norm(B_eci(:,1));

%Normalize the satellite measurements
B_meas(:,1) = B_sat(1,:)'/ norm(B_eci(:,1));
W_meas(:,1) = qw(5:7);
 
%Loop responsible for the integration of the attitude motion.
%The estimated quaternion is used to rotate the magnectic induction vector
%provided by the IGRF model for each time step

for i = 2 : 1 : length(t)
    %Notice that the attitude motion integration here is used for the cost
    %function determination, not for the actual kalman filter state
    %estimation

% 	A = @(q)((q(1)^2 - norm(q(2:4))^2)*eye(3) - 2*q(1)*[0, -q(4), q(3); q(4), 0, -q(2); -q(3), q(2), 0] + 2*q(2:4)*q(2:4)');
    A = @quat2dcm;

% State equations
% first line = attitude kinematics = 1/2*q o w
% second line = attitude dynamics = Jdw + w x Jw = Mgg
	f = @(t,x) [-0.5*x(5:7)'*x(2:4);
                0.5*(x(5)*x(1)+x(7)*x(3)-x(6)*x(4));
                0.5*(x(6)*x(1)-x(7)*x(2)+x(5)*x(4));
                0.5*(x(7)*x(1)+x(6)*x(2)-x(5)*x(3));
%         (0.5*[x(4), -x(3), x(2); x(3), x(4), -x(1); -x(2), x(1), x(4); -x(1), -x(2), -x(3)]*x(5:7));
            J\(3*(worb^2)*cross(A(x(1:4)')*e3, J*A(x(1:4)')*e3) -cross(x(5:7), J*x(5:7)) +cross(x(8:10),B_model(:, i-1))*1e-9 );
            zeros(3,1);
            zeros(3,1)];


% Numerical integration
    [~, res] = ode45(f,[t(i-1) t(i)], qw); 
    
    qw = res(end,:)';
    e3 = xsat_eci(:,i)/norm(xsat_eci(:,i));
    W_meas(:, i) = qw(5:7); % + qw(11:13);
    
    %Determine the IGRF model for the sattelite attitude in the respective
    %sample of data
    
    B_model(:, i) = quatrotate(qw(1:4)', B_eci(:, i)')' + qw(11:13);    
    
    %Normalization of both IGRF model and IMU measurements model 
    B_model_norm(:, i) = B_model(:, i)/ norm(B_eci(:, i));             
    B_meas(:, i) = B_sat(i, :)'/ norm(B_eci(:, i));

end


% Cost function computation
F = sum((B_meas(1,:)-B_model_norm(1,:)).^2+(B_meas(2,:)-B_model_norm(2,:)).^2+(B_meas(3,:)-B_model_norm(3,:)).^2)
% F = F+sum((W_meas(1,:)-Omega_meas(1,:)).^2+(W_meas(2,:)-Omega_meas(2,:)).^2+(W_meas(3,:)-Omega_meas(3,:)).^2); 

%-----------------------------------------------------------------------------
%% Plots of the Differential Evolution results 
%-----------------------------------------------------------------------------


figure('Color',[1 1 1])
plot(t,B_meas(1,:),'b','LineWidth',2)
hold on
plot(t,B_meas(2,:),'g','LineWidth',2)
plot(t,B_meas(3,:),'r','LineWidth',2)
plot(t,B_model_norm(1,:),'--b','LineWidth',2)
plot(t,B_model_norm(2,:),'--g','LineWidth',2)
plot(t,B_model_norm(3,:),'--r','LineWidth',2)
grid on
xlabel('Tempo [s]')
ylabel('Componentes do vetor campo magnético unitário')
legend('B_x measurements','B_y measurements','B_z measurements','B_x model','B_y model','B_z model')



figure('Color',[1 1 1])
plot(t,Omega_meas(1,:)*180/pi,'b','LineWidth',2)
hold on
plot(t,Omega_meas(2,:)*180/pi,'g','LineWidth',2)
plot(t,Omega_meas(3,:)*180/pi,'r','LineWidth',2)
plot(t,W_meas(1,:)*180/pi,'--b','LineWidth',2)
plot(t,W_meas(2,:)*180/pi,'--g','LineWidth',2)
plot(t,W_meas(3,:)*180/pi,'--r','LineWidth',2)
grid on
xlabel('Tempo [s]')
ylabel('Velocidade angular [graus/s]')
legend('\omega_x meas','\omega_y meas','\omega_z meas','\omega_x model','\omega_y model','\omega_z model')

W_shift = [-1.4;-2.9; -0.5]; % July 9
% W_shift = [-1.4;-1.9; -0.5]; % August 7
figure('Color',[1 1 1])
plot(t,Omega_meas(1,:)*180/pi-W_shift(1),'b','LineWidth',2)
hold on
plot(t,Omega_meas(2,:)*180/pi-W_shift(2),'g','LineWidth',2)
plot(t,Omega_meas(3,:)*180/pi-W_shift(3),'r','LineWidth',2)
plot(t,W_meas(1,:)*180/pi,'--b','LineWidth',2)
plot(t,W_meas(2,:)*180/pi,'--g','LineWidth',2)
plot(t,W_meas(3,:)*180/pi,'--r','LineWidth',2)
grid on

xlabel('Time, s')
ylabel('Angular velocity, deg/s')
legend('\omega_x measurements without shift','\omega_y measurements without shift','\omega_z measurements without shift','\omega_x model','\omega_y model','\omega_z model')



end


