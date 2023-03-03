 function [omega_orbital]= orbangvel( xsat_eci, vsat_eci);
  
omega_eci=cross(xsat_eci, vsat_eci)/((xsat_eci(1))^2+(xsat_eci(2))^2+(xsat_eci(3))^2);


Z_orb=xsat_eci/sqrt(xsat_eci(1)^2+xsat_eci(2)^2+xsat_eci(3)^2);
X1_orb=vsat_eci/sqrt(vsat_eci(1)^2+vsat_eci(2)^2+vsat_eci(3)^2);
Y_orb=cross(Z_orb, X1_orb);
Y_orb=Y_orb/sqrt(Y_orb(1)^2+Y_orb(2)^2+Y_orb(3)^2);
X_orb=cross(Y_orb, Z_orb);
X_orb=X_orb/sqrt(X_orb(1)^2+X_orb(2)^2+X_orb(3)^2);


omega_orbital(1,1)=omega_eci(1)*X_orb(1)+omega_eci(2)*X_orb(2)+omega_eci(3)*X_orb(3);
omega_orbital(2,1)=omega_eci(1)*Y_orb(1)+omega_eci(2)*Y_orb(2)+omega_eci(3)*Y_orb(3);
omega_orbital(3,1)=omega_eci(1)*Z_orb(1)+omega_eci(2)*Z_orb(2)+omega_eci(3)*Z_orb(3);
 end
