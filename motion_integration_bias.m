function [A_quat, qw] = motion_integration_bias(qw, x_eci, B_model, index)
% Motion Integration for the 

global t J worb

% Local Vertical Vector in ECI frame
e3 = x_eci/norm(x_eci);

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
            J\(3*(worb^2)*cross(A(x(1:4)')*e3, J*A(x(1:4)')*e3) -cross(x(5:7), J*x(5:7)) + cross(x(8:10),B_model)*1e-9);
            zeros(3,1);
            zeros(3,1)];


% numerical integration
    [~, res] = ode45(f,[t(index) t(index + 1)], qw); 
    
    qw = res(end,:)';
    
    A_quat = A(qw(1:4)');
    
end


