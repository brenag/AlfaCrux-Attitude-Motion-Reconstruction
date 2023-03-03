function [F] = cost_func(qw)


global B_eci B_sat t J worb xsat_eci
%Notice that the 'for loop' responsible for the integration of the attitude
%motion starts in index 2, so we must determine the magnetic field
%reference for the initial time
e3 = xsat_eci(:,1)/norm(xsat_eci(:,1));
%Determine the IGRF model for the sattelite initial quaternion guess
B_model(:, 1) = quatrotate(qw(1:4)', B_eci(:,1)')'+qw(11:13);
B_model_norm(:, 1) = B_model(:, 1)/ norm(B_eci(:,1));

%Normalize the satellite measurements
B_meas(:,1) = B_sat(1,:)'/ norm(B_eci(:,1)); 

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
% numerical integration
    
% numerical integration
    [~, res] = ode45(f,[t(i-1) t(i)], qw);
    
    qw = res(end,:)';
    
    %Determine the IGRF model for the sattelite attitude in the respective
    %sample of data
    e3 = xsat_eci(:,i)/norm(xsat_eci(:,i));
    B_model(:, i) = quatrotate(qw(1:4)', B_eci(:, i)')'+qw(11:13);    
    
    %Normalization of both IGRF model and IMU measurements model 
    B_model_norm(:, i) = B_model(:, i)/ norm(B_eci(:, i));             
    B_meas(:, i) = B_sat(i, :)'/ norm(B_eci(:, i));

end

%Loop responsible to calculate the cost function using the magnectic
%induction vector measured by the IMU and estimated by the IGRF model

%Alternative expression for the cost function without the 'for loop'
F = sum((B_meas(1,:)-B_model_norm(1,:)).^2+(B_meas(2,:)-B_model_norm(2,:)).^2+(B_meas(3,:)-B_model_norm(3,:)).^2); 
end


