 function [S_orbital, S_eci]= Sunvec(sat_llh, jdnow, xsat_eci, vsat_eci);
  
 [year_now,mon_now,day_now,hr_now,min_now,sec_now,days] = invjday ( jdnow );
 
 
    nd=days;
    UT=hr_now+min_now/60+sec_now/3600;
    U0=12;
    K=15;
    alpha2=23.5*pi/180;
    
    phiSE=360*(172-nd)/365*pi/180;
    beta=asin(sin(alpha2)*cos(phiSE));
   
    
    
    %%солнечный ветер
%     r1=10000*(np+4*na)^(-1/6)*V^(-1/3)*Rearth;
    
    
beta1=K*(UT-U0)*pi/180;

   
       
phideg=sat_llh(1)*180/pi;
h=sat_llh(3);
lambda=sat_llh(2);

phi=phideg*pi/180;

   %%учёт эллиптичности Земли

a=6378200;
b=6356800;
Rearth=6371200;  

R=sqrt((h^2)+2*h*sqrt(((a^2)*(cos(phi))^2)+(b^2)*(sin(phi))^2)+(((a^4)*(cos(phi))^2)+(b^4)*((sin(phi))^2))/(((a^2)*(cos(phi))^2)+(b^2)*(sin(phi))^2));
phi1=atan((((b^2)+h*sqrt(((a^2)*(cos(phi))^2)+(b^2)*(sin(phi))^2))/(((a)^2)+h*sqrt(((a^2)*(cos(phi))^2)+(b^2)*(sin(phi))^2))*tan(phi)));
teta=pi/2-phi1;



    
%%матрица перехода из геоцентрических к солнечно-магнитосферным коордиранатам
    
% T=[cos(beta1)*cos(beta)  -sin(beta1)*cos(beta)  sin(beta); sin(beta1)*cos(beta2)-cos(beta1)*sin(beta)*sin(beta2)  cos(beta1)*cos(beta2)+sin(beta1)*sin(beta)*sin(beta2)  cos(beta)*sin(beta2); -sin(beta1)*sin(beta2)-cos(beta1)*sin(beta)*cos(beta2)  -cos(beta1)*sin(beta2)+sin(beta1)*sin(beta)*cos(beta2)  cos(beta)*cos(beta2)];
    

%%матрица перехода из сферических к экваториальным координатам
   
S=[sin(teta)*cos(lambda)  cos(teta)*cos(lambda)  -sin(lambda); sin(teta)*sin(lambda)  cos(teta)*sin(lambda)  cos(lambda); cos(teta)  -sin(teta)  0];
% S=[cos(phi1)*cos(lambda)  sin(phi1)*cos(lambda)  -sin(lambda); cos(phi1)*sin(lambda)  sin(phi1)*sin(lambda)  cos(lambda); sin(phi1)  -cos(phi1)  0];

% SolMagn=T*Equat;

%% Вектор направления на Солнце в Гринвичских координатах 

S_equat(1,1)=cos(beta1)*cos(beta);
S_equat(2,1)=(-1)*sin(beta1)*cos(beta);
S_equat(3,1)=sin(beta);
        
%% Приведение к одной СК и вывод компонент вектора направления на Солнце в орбитальной СК



S1=gstime(jdnow);
P=[   cos(S1) -sin(S1) 0;
      sin(S1) cos(S1)  0;
        0        0     1]; 
S_geocentr=P*S_equat;
S_eci=S_geocentr;

Z_orb=xsat_eci/sqrt(xsat_eci(1)^2+xsat_eci(2)^2+xsat_eci(3)^2);
X1_orb=vsat_eci/sqrt(vsat_eci(1)^2+vsat_eci(2)^2+vsat_eci(3)^2);
Y_orb=cross(Z_orb, X1_orb);
Y_orb=Y_orb/sqrt(Y_orb(1)^2+Y_orb(2)^2+Y_orb(3)^2);
X_orb=cross(Y_orb, Z_orb);
X_orb=X_orb/sqrt(X_orb(1)^2+X_orb(2)^2+X_orb(3)^2);


S_orbital(1,1)=S_geocentr(1)*X_orb(1)+S_geocentr(2)*X_orb(2)+S_geocentr(3)*X_orb(3);
S_orbital(2,1)=S_geocentr(1)*Y_orb(1)+S_geocentr(2)*Y_orb(2)+S_geocentr(3)*Y_orb(3);
S_orbital(3,1)=S_geocentr(1)*Z_orb(1)+S_geocentr(2)*Z_orb(2)+S_geocentr(3)*Z_orb(3);
 end
