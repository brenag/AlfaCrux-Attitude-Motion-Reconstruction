% Входные данные: 
% sat_llh    - ???
% jdnow      - время в эпохе j2000
% xsat_eci   - координаты спутника в ИСК
% vsat_eci   - вектор скорости спутника в ИСК
%
% Выходные данные:
% S_orbital  - направление на Солнце в ОСК
% S_eci      - направление на Солнце в ИСК
% 
%  See also Sunvec, embed_func, ReadTLE, ReadTLE_name
%  
function [S_orbital, S_eci]= Sunvec_1(jdnow, xsat_eci, vsat_eci)
  
    %%  
    %     Values determined using data from 1950-1991 in the 1990 Astronomical
    %     Almanac.  See DELTA_ET.WQ1 for details. 
    
    %S_eci - координаты направления на Солнце в ИСК (если не так, то ПЕРЕКОДИТЬ ВСЁ;ь)
    %S_orbital - соответственно, в ОСК
    
    mjd  = jdnow - 2415020.0;
    year = 1900 + mjd/365.25;
%     SECDAY=rem(jdnow,1)*3600;
    SECDAY=86400;
    T = (mjd + Delta_ET(year)/SECDAY)/36525.0;
    M = Modulus(358.47583 + Modulus(35999.04975*T,360.0)-... 
      (0.000150 + 0.0000033*T)*(T^2),360.0) * pi/180;
    L = Modulus(279.69668 + Modulus(36000.76892*T,360.0)+... 
      0.0003025*(T^2),360.0) * pi/180;
    e = 0.01675104 - (0.0000418 + 0.000000126*T)*T;
    C = ((1.919460 - (0.004789 + 0.000014*T)*T)*sin(M)+... 
       (0.020094 - 0.000100*T)*sin(2*M) + 0.000293*sin(3*M))*pi/180;
    O = (Modulus(259.18 - 1934.142*T,360.0))*pi/180;
    Lsa = Modulus(L + C - pi/180*(0.00569 - 0.00479*sin(O)),2*pi);
    nu  = Modulus(M + C,2*pi);
    R  = 1.0000002*(1 - (e^2))/(1 + e*cos(nu));
    eps = (23.452294 - (0.0130125 + (0.00000164 - 0.000000503*T)*T)*T+... 
                0.00256*cos(O))*pi/180;
    %  R  = AU*R;
    S_eci(1) = R*cos(Lsa);
    S_eci(2) = R*sin(Lsa)*cos(eps);
    S_eci(3) = R*sin(Lsa)*sin(eps);
    
    %тут дальше как было в Sunvec
    
    Z_orb=xsat_eci/sqrt(xsat_eci(1)^2+xsat_eci(2)^2+xsat_eci(3)^2);
    X1_orb=vsat_eci/sqrt(vsat_eci(1)^2+vsat_eci(2)^2+vsat_eci(3)^2);
    Y_orb=cross(Z_orb, X1_orb);
    Y_orb=Y_orb/sqrt(Y_orb(1)^2+Y_orb(2)^2+Y_orb(3)^2);
    X_orb=cross(Y_orb, Z_orb);
    X_orb=X_orb/sqrt(X_orb(1)^2+X_orb(2)^2+X_orb(3)^2);


    S_orbital(1)=S_eci(1)*X_orb(1)+S_eci(2)*X_orb(2)+S_eci(3)*X_orb(3);
    S_orbital(2)=S_eci(1)*Y_orb(1)+S_eci(2)*Y_orb(2)+S_eci(3)*Y_orb(3);
    S_orbital(3)=S_eci(1)*Z_orb(1)+S_eci(2)*Z_orb(2)+S_eci(3)*Z_orb(3);
    
end

% void CSunModel::CalcSolarPos( double JulianDate,  vector& SunPosition )
% {
%   double mjd,year,T,M,L,e,C,O,Lsa,nu,R,eps;
% //  double ob[5];
%   mjd  = JulianDate - 2415020.0;
%   year = 1900 + mjd/365.25;
%   T = (mjd + Delta_ET(year)/SECDAY)/36525.0;
%   M = Modulus(358.47583 + Modulus(35999.04975*T,360.0)- 
% 	  (0.000150 + 0.0000033*T)*pow(T,2),360.0) * pi/180;
%   L = Modulus(279.69668 + Modulus(36000.76892*T,360.0)+ 
% 	  0.0003025*pow(T,2),360.0) * pi/180;
%   e = 0.01675104 - (0.0000418 + 0.000000126*T)*T;
%   C = ((1.919460 - (0.004789 + 0.000014*T)*T)*sin(M)+ 
% 	   (0.020094 - 0.000100*T)*sin(2*M) + 0.000293*sin(3*M))*pi/180;
%   O = (Modulus(259.18 - 1934.142*T,360.0))*pi/180;
%   Lsa = Modulus(L + C - pi/180*(0.00569 - 0.00479*sin(O)),2*pi);
%   nu  = Modulus(M + C,2*pi);
%   R  = 1.0000002*(1 - pow(e,2))/(1 + e*cos(nu));
%   eps = (23.452294 - (0.0130125 + (0.00000164 - 0.000000503*T)*T)*T+ 
% 				0.00256*cos(O))*pi/180;
% //  R  = AU*R;
%   SunPosition.x = R*cos(Lsa);
%   SunPosition.y = R*sin(Lsa)*cos(eps);
%   SunPosition.z = R*sin(Lsa)*sin(eps);
% }

%% Modulus

%?????????????????????????????????????????????????????????
%?????????????????????????????????????????????????????????
%???????modu = arg1 - (long)(arg1/arg2) * arg2;???????????
%?????????????????????????????????????????????????????????
%?????????????????????????????????????????????????????????

function [M]= Modulus(arg1, arg2)
%double CSunModel::Modulus(double arg1, double arg2)
    %???????????????????????????????????????
    %modu = arg1 - (long)(arg1/arg2) * arg2;
    %fix - округление до ближайшего целого
    modu = arg1 - fix(arg1/arg2) * arg2;
	if modu >= 0	
        M = modu;
    else M = modu + arg2;
    end
end

% double CSunModel::Modulus(double arg1, double arg2)
% {
% 	double modu;
% 	modu = arg1 - (long)(arg1/arg2) * arg2;
% 	if (modu >= 0.0)	return modu;
% 	else			return modu + arg2;
% }

%% Delta_ET
function [a]= Delta_ET(year);
% double CSunModel::Delta_ET( double year)
    %  /* { Values determined using data from 1950-1991 in the 1990 Astronomical
    %     Almanac.  See DELTA_ET.WQ1 for details. }*/
    a = (26.465 + 0.747622*(year - 1950) + 1.886913*sin(2*pi*(year - 1975)/33));
    %//{Function Delta_ET}
end

% double CSunModel::Delta_ET( double year)
% {
%  /* { Values determined using data from 1950-1991 in the 1990 Astronomical
%     Almanac.  See DELTA_ET.WQ1 for details. }*/
%   return (26.465 + 0.747622*(year - 1950) + 1.886913*sin(2*pi*(year - 1975)/33));
% };
% //{Function Delta_ET}

