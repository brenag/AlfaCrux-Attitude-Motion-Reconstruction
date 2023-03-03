function [F] = Func(B_bias)
% Function to be minimized by fminunc
    global B

    % Mod variable is the average of the modulus of the IGRF magnetic field
    % during the period of the analysis

    Mod=20750;
    
    Bx = B_bias(1);
    By = B_bias(2);
    Bz = B_bias(3);
    Mod_exp = ((B(:,1)-Bx).^2 + (B(:,2)-By).^2 + (B(:,3)-Bz).^2).^(1/2);

    Diff=Mod_exp-Mod;
    F=sum(Diff.^2);

end

