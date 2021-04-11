% ---------------------------------- %
%  Potential transpiration partition %
% -----------------------------------%
function [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, D50, c, b, theta_sat, wet, zm)
    % ------- function input -------
    % pEc    :  Potnetial Evaporation on canopy
    % w      :  Initialized values for soil moisture
    % PFTpar :  PFT parameters
    % ------- function output ------
    % pTr_ly :  seperate potnetial Transpiration
    % -------

    z1 = zm(1);
    z2 = zm(2);
    z3 = zm(3);
    
    r1 = (1 / (1 + (z1 / D50)^c)); 
    r2 = (1 / (1 + (z2 / D50)^c)) - (1 / (1 + (z1 / D50)^c)); 
    r3 = (1 / (1 + (z3 / D50)^c)) - (1 / (1 + (z2 / D50)^c)); 
    
    % the maximum transpiration rate of each soil layer, Tr_p
    % Get all avaliable water contents through root distribution
    wr = r1 * (wa1 / theta_sat)^b + r2 * (wa2 / theta_sat)^b + r3 * (wa3 / theta_sat)^b;

    % Root distribution adjusted by soil water content
    beta1 = r1 * (wa1 / theta_sat)^b / wr;
    beta2 = r2 * (wa2 / theta_sat)^b / wr;
    beta3 = r3 * (wa3 / theta_sat)^b / wr;

    % potentail transpiration rate for different layers
    Tr_p1 = (1 - wet) * beta1 * pEc;
    Tr_p2 = (1 - wet) * beta2 * pEc;
    Tr_p3 = (1 - wet) * beta3 * pEc;

end