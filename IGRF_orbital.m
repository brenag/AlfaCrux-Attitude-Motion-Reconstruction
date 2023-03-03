 function [B_orbital, B_eci]= IGRF_orbital(sat_llh, jdnow, day_now, hr_now, min_now, sec_now, xsat_eci, vsat_eci, n1);
 % Secular Variation

load('IGRF_coeffs.mat');

jd_2010 = jday(2010,  1,  1,  0,  0,  0);
jd_2015 = jday(2014, 12, 31, 23, 59, 59);

Gnm=Gnm_vvod+deltaGnm*((jdnow-jd_2010)/(jd_2015-2010)*5);
Hnm=Hnm_vvod+deltaHnm*((jdnow-jd_2010)/(jd_2015-2010)*5);

%angles
    [year_now,mon_now,day_now,hr_now,min_now,sec_now,days] = invjday ( jdnow );
 
 
    nd=days;
    UT=hr_now+min_now/60+sec_now/3600;
    U0=12;
    K=15;
    alpha2=23.5*pi/180;
    alpha1=11*pi/180;
    
    phiSE=360*(172-nd)/365*pi/180;
    phim=(K*UT-69)*pi/180;
    beta=asin(sin(alpha2)*cos(phiSE));
    psi=asin(-sin(beta)*cos(alpha1)+cos(beta)*sin(alpha1)*cos(phim));
    
    
    %%солнечный ветер
%     r1=10000*(np+4*na)^(-1/6)*V^(-1/3)*Rearth;
    
    
beta1=K*(UT-U0)*pi/180;
beta2=acos((cos(alpha1)+sin(beta)*sin(psi))/(cos(psi)*cos(beta))); 
   
       
phideg=sat_llh(1)*180/pi;
h=sat_llh(3);
lambda=sat_llh(2);


%%широта на экваторе и у полюсов
if phideg>89.99999999999999289 & phideg<=90
    phideg=89.99999999999999289;   end
if phideg>=0 & phideg<(1.2164*10^-13)
    phideg=1.2164*10^-13;          end
if phideg<0 & phideg>(-1.2163*10^-13)
    phideg=-1.2163*10^-13;         end


phi=phideg*pi/180;

   %%Elliptial model of Earth

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



%% Магнитное поле Земли внутренних источников


Bx=0;
By=0;
Bz=0;


for n=1:n1
for m=0:n
    
    m1=m+1; %%небольшое изменение нумерации для матриц Gnm и Hnm
    
    if m<1  epsm=1;
    else    epsm=2;
    end
    
    summ1=0;
    summ2=0;
    
    mnozhitel1=1;
    mnozhitel2=1;
    mnozhitel3=1;  
    
    
    for s=1:n
    mnozhitel1=mnozhitel1*(2*s-1);
    end;
    
    for k=1:n1
        
        mnozhitel2=mnozhitel2*(n-m+2-2*k)*(n-m+1-2*k);
        mnozhitel3=mnozhitel3*(2*k*(2*n+1-2*k));
    
        
        summ1=summ1+((-1)^k)*(mnozhitel2/mnozhitel3)*(cos(teta))^(n-m-2*k);
        summ2=summ2+((-1)^k)*(mnozhitel2/mnozhitel3)*(n-m-2*k)*(cos(teta))^(n-m-2*k-1);
        
    end;
    
    
    Bx=Bx+(Gnm(n,m1)*cos(m*lambda)+Hnm(n,m1)*sin(m*lambda))*(mnozhitel1*sqrt(epsm/(factorial(n+m)*factorial(n-m))))*((sin(teta))^m)*((m*(cos(teta)/sin(teta))*((cos(teta))^(n-m)+summ1))-sin(teta)*((n-m)*(cos(teta))^(n-m-1)+summ2))*(Rearth/R)^(n+2); 
    By=By+m*(Gnm(n,m1)*sin(m*lambda)-Hnm(n,m1)*cos(m*lambda))*(mnozhitel1*sqrt(epsm/(factorial(n+m)*factorial(n-m)))*((sin(teta))^m)*((cos(teta))^(n-m)+summ1))/sin(teta)*(Rearth/R)^(n+2);
    Bz=Bz-(n+1)*(Gnm(n,m1)*cos(m*lambda)+Hnm(n,m1)*sin(m*lambda))*(mnozhitel1*sqrt(epsm/(factorial(n+m)*factorial(n-m)))*((sin(teta))^m)*((cos(teta))^(n-m)+summ1))*(Rearth/R)^(n+2);
    
    
end;    
end;

    
%% Магнитосфера Земли 
        
%     
%     
%     F=[Q(1,1)*sin(psi); Q(2,1)*sin(psi)*cos(psi); Q(3,1)*sin(psi); Q(4,1)*(cos(psi))^2+Q(5,1)*(sin(psi))^2; Q(6,1)*cos(psi); Q(7,1)*(cos(psi))^2+Q(8,1)*(sin(psi))^2; Q(9,1)*cos(psi); Q(10,1)*sin(psi)*cos(psi)];
%     G=[S(1,1)*cos(psi); S(2,1); S(1,1)*sin(psi)];
%     H=[-Q(1,1)*cos(psi); -Q(4,1)*(sin(psi))^2-Q(5,1)*(cos(psi))^2; -Q(3,1)*cos(psi); -Q(2,1)*sin(psi)*cos(psi); Q(6,1)*sin(psi); Q(10,1)*sin(psi)*cos(psi); Q(9,1)*sin(psi); Q(7,1)*(sin(psi))^2+Q(8,1)*(cos(psi))^2];
%     
%     B2x=F(1,1)+F(2,1)*SolMagn(1,1)/r1+F(3,1)*SolMagn(2,1)/r1+F(4,1)*SolMagn(3,1)/r1+((psi*180/pi)/10)*(F(5,1)+F(6,1)*SolMagn(1,1)/r1+F(7,1)*SolMagn(2,1)/r1+F(8,1)*SolMagn(3,1)/r1);
%     B2y=((psi*180/pi)/10)*(G(1,1)*SolMagn(1,1)/r1+G(2,1)*SolMagn(2,1)/r1+G(3,1)*SolMagn(3,1)/r1);
%     B2z=H(1,1)+H(2,1)*SolMagn(1,1)/r1+H(3,1)*SolMagn(2,1)/r1+H(4,1)*SolMagn(3,1)/r1+((psi*180/pi)/10)*(H(5,1)+H(6,1)*SolMagn(1,1)/r1+H(7,1)*SolMagn(2,1)/r1+H(8,1)*SolMagn(3,1)/r1);
    
 
    
    
    
%% Приведение к одной СК и вывод компонент вектора напряжённости магнитного поля в орбитальной СК


B1GeoSphere=[-Bz; -Bx; By];

B_equat=S*B1GeoSphere;


S1=gstime(jdnow);
P=[   cos(S1) -sin(S1) 0;
      sin(S1) cos(S1)  0;
        0        0     1]; 
B_geocentr=P*B_equat;
B_eci=B_geocentr;


Z_orb=xsat_eci/sqrt(xsat_eci(1)^2+xsat_eci(2)^2+xsat_eci(3)^2);
X1_orb=vsat_eci/sqrt(vsat_eci(1)^2+vsat_eci(2)^2+vsat_eci(3)^2);
Y_orb=cross(Z_orb, X1_orb);
Y_orb=Y_orb/sqrt(Y_orb(1)^2+Y_orb(2)^2+Y_orb(3)^2);
X_orb=cross(Y_orb, Z_orb);
X_orb=X_orb/sqrt(X_orb(1)^2+X_orb(2)^2+X_orb(3)^2);


B_orbital(1,1)=B_geocentr(1)*X_orb(1)+B_geocentr(2)*X_orb(2)+B_geocentr(3)*X_orb(3);
B_orbital(2,1)=B_geocentr(1)*Y_orb(1)+B_geocentr(2)*Y_orb(2)+B_geocentr(3)*Y_orb(3);
B_orbital(3,1)=B_geocentr(1)*Z_orb(1)+B_geocentr(2)*Z_orb(2)+B_geocentr(3)*Z_orb(3);

 end
