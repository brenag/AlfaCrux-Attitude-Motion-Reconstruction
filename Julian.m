function [JD] = Julian(yr,mo,d,h,min,s)
A=double(2-int32(yr/100)+int32(int32(yr/100)/4));
B=double(int32(365.25*(yr+4716))+int32(30.6001*(mo+1)));
JD=B+A+d-1524.5+(((s/60+min)/60+h)/24);
end

