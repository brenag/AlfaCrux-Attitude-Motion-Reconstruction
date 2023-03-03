function [omega_orbital, S_orbital, B_orbital, S_eci, B_eci, xsat_eci, vsat_eci, gst, SHADOW, sat_llh, cos_delta, a, xsat_ecf, vsat_ecf]= embed_func(t, t_begin, flag1, flag2, flag3, longstr1, longstr2);
% Function with SGP4, IGRF, and Sun direction models


% Extract the Orbital Elements from TLE to use in SGP4 
satrec = twoline2rvMOD(longstr1,longstr2);

c=t_begin;
startyear = c(1);
startmon  = c(2);
startday  = c(3);
starthr   = c(4);
startmin  = c(5);
startsec  = c(6);

% Gets Epoch (time of the first sample) in Julian Date
jdstart = jday( startyear,startmon,startday,starthr,startmin,startsec );

[ startmfe ] = timing(jdstart, satrec); % time past since TLE Epoch until start time

[ jdnow  ] = timing2(jdstart, t); % Julian Date of current sample

[year_now,mon_now,day_now,hr_now,min_now,sec_now] = invjday (jdnow); 

[ nowmfe ] = timing(jdnow, satrec); % time past since TLE Epoch until current sample

% Initiate some variables 
gst = 0 ;
xsat_eci = [0; 0; 0];
vsat_eci = [0; 0; 0]; 

% SGP4 propagation in ECI and conversion to ECEF to use in IGRF
[satrec, xsat_ecf, vsat_ecf, xsat_eci, vsat_eci, gst] = sgp4_ecf(satrec, nowmfe);

   xsat_ecf=xsat_ecf*1000;  %m
   vsat_ecf=vsat_ecf*1000;  %mps
   
   xsat_eci=xsat_eci*1000;  %m
   vsat_eci=vsat_eci*1000;  %mps

   % Obtain lattitude, longitude and altitude of current position
   sat_llh = [0; 0; 0];
   sat_llh = ecf2llhT(xsat_ecf); 

% Obtain the Geomagnetic induction vector using IGRF model
if flag2==1
    n1=13;
    [B_orbital, B_eci] = IGRF_orbital (sat_llh, jdnow, day_now, hr_now, min_now, sec_now, xsat_eci, vsat_eci, n1);
end

% Obtain the Sun direction vector
if flag3==1
[S_orbital, S_eci] = Sunvec (sat_llh, jdnow, xsat_eci, vsat_eci);
elseif flag3==2
[S_orbital, S_eci] = Sunvec_1(jdnow, xsat_eci, vsat_eci);
end


% Some calculations to find out if the satellite is on eclipse
cos_delta=(xsat_eci(1)*S_eci(1)+xsat_eci(2)*S_eci(2)+xsat_eci(3)*S_eci(3))/sqrt(xsat_eci(1)^2+xsat_eci(2)^2+xsat_eci(3)^2);
a=sqrt(xsat_eci(1)^2+xsat_eci(2)^2+xsat_eci(3)^2);

if flag1==1
    if cos_delta<0 
        delta=acos(cos_delta);
        if  abs(sqrt(xsat_eci(1)^2+xsat_eci(2)^2+xsat_eci(3)^2)*cos(delta-pi/2)/1000) <= 6371.2
        SHADOW=1;
        else
        SHADOW=0;    
        end
   
    else
        SHADOW=0;      
    end
else 
SHADOW=0;
end

% Orbital Velocity
[omega_orbital]= orbangvel(xsat_eci, vsat_eci);


end

