function [mfe  ] = timing(jd, satrec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

            mfe = (jd - satrec.jdsatepoch) * 1440;

% %    
end

