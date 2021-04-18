% --------------- %
%      Runoff     %
% ----------------%
function [srf, IWS] = runoff(Pnet, zgw, zm, wa, soilpar)
    % --------- function input -------
    % zm       : Soil depth of different layers
    % wa       : The antecedent soil water content expressed
    %            as a function of the WHC in that layer
    % soilpar  : Soil parameters according to Soil type
    % Pnet     : Net precipitation = P-I+Snowmelt
    % zgw      : groundwater table depth
    % --------- function output ------
    % srf      : Surface Runoff, mm
    % IWS      : Water enter into soil surface, mm
    % ---------
    % Reference:
    % Choudhury BJ, Digirolamo NE, 1998
    % SCS (1985). National Engineering Handbook. Section 4: Hydrology.
    % Washington, DC: Soil Conservation Service, U.S. Department of Agriculture.
    % SCS (1986). Urban Hydrology for Small Watersheds, Technical Release No. 55.
    % Washington, DC: Soil Conservation Service, U.S. Department of Agriculture.
    % -------------------------------------------------------------------------

    % saturated wa for specific soil type
    theta_sat = soilpar(3);

    if zgw <= 0

        % exceeded groundwater on the soil surface
        Vmax = 0;
        srf = 0 - zgw * theta_sat;

    elseif zgw > 0 && zgw <= zm(1)

        % the thickness of unsaturated soil in ith layer, (mm)
        d1 = zgw;

        % the unsaturated soil water in ith layer, (mm)
        wa1_unsat = (wa(1) * zm(1) - theta_sat * (zm(1) - d1)) / d1;

        % calculate the overall soil water retention capacity, Vmax
        Vmax = d1 * (theta_sat - wa1_unsat);

        if Pnet > 0.2 * Vmax
            srf = (Pnet - 0.2 * Vmax)^2 / (Pnet + 0.8 * Vmax);
        else
            srf = 0;
        end

    elseif zgw > zm(1) && zgw <= zm(1) + zm(2)

        % the thickness of unsaturated soil in ith layer, (mm)
        d1 = zm(1);
        d2 = zgw - zm(1);

        % the unsaturated soil water in ith layer, (mm)
        wa1_unsat = wa(1);
        wa2_unsat = (wa(2) * zm(2) - theta_sat * (zm(2) - d2)) / d2;

        % calculate the overall soil water retention capacity, Vmax
        Vmax = d1 * (theta_sat - wa1_unsat) + d2 * (theta_sat - wa2_unsat);

        if Pnet > 0.2 * Vmax
            srf = (Pnet - 0.2 * Vmax)^2 / (Pnet + 0.8 * Vmax);
        else
            srf = 0;
        end

    elseif zgw > zm(1) + zm(2) && zgw <= zm(1) + zm(2) + zm(3)

        % the thickness of unsaturated soil in ith layer, (mm)
        d1 = zm(1);
        d2 = zm(2);
        d3 = zgw - zm(2) - zm(1);

        % the unsaturated soil water in ith layer, (mm)
        wa1_unsat = wa(1);
        wa2_unsat = wa(2);
        wa3_unsat = (wa(3) * zm(3) - theta_sat * (zm(3) - d3)) / d3;

        % calculate the overall soil water retention capacity, Vmax
        Vmax = d1 * (theta_sat - wa1_unsat) + d2 * (theta_sat - wa2_unsat) + ...
            d3 * (theta_sat - wa3_unsat);

        if Pnet > 0.2 * Vmax
            srf = (Pnet - 0.2 * Vmax)^2 / (Pnet + 0.8 * Vmax);
        else
            srf = 0;
        end

    elseif zgw > zm(1) + zm(2) + zm(3)

        % the thickness of unsaturated soil in ith layer, (mm)
        d1 = zm(1);
        d2 = zm(2);
        d3 = zm(3);

        % the unsaturated soil water in ith layer, (mm)
        wa1_unsat = wa(1);
        wa2_unsat = wa(2);
        wa3_unsat = wa(3);

        % calculate the overall soil water retention capacity, Vmax
        Vmax = d1 * (theta_sat - wa1_unsat) + d2 * (theta_sat - wa2_unsat) + ...
            d3 * (theta_sat - wa3_unsat);

        if Pnet > 0.2 * Vmax
            srf = (Pnet - 0.2 * Vmax)^2 / (Pnet + 0.8 * Vmax);
        else
            srf = 0;
        end
        
    end

    % actual water enter into the soil surface, (mm)
    IWS = min(Vmax, Pnet - srf);
    IWS = max(0, IWS);

    % redundant water --> runoff (balance)
    if IWS == Vmax || Vmax <= 0
        srf = Pnet - Vmax;
    end

end
