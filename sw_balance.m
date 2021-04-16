% ----------------------- %
%    Soil Water Balance   %
% ----------------------- %
function [wa, zgw, Tr, Es, uex] = sw_balance(IWS, pEc, pEs, Ta, Topt, wa, soilpar, pftpar, wet, zm, zgw)
    % ---------- function input -------
    % IWS     -- total water enter into soil surface, mm
    % pEc     -- potential ET allocate to plant, mm
    % pEs     -- potential ET allocate to soil surface, mm
    % Ta      -- air temperature, C
    % Topt    -- optimal growth temperature for plant, C
    % wa      -- previous soil water content, 3 layers
    % soilpar -- soil-related parameters
    % pftpar  -- plant-related parameters
    % wet     -- wetness fraction indice
    % zm      -- soil layer depth, 3 layers
    % zgw     -- groundwater table depth, mm 
    % ---------- function output -------
    % Tr      -- actual plant transpiration, mm
    % Es      -- actual soil evaporation, mm
    % wa      -- updated soil water content, 3 layers, %
    % zgw     -- groundwater table depth, mm
    % uex     -- goundwater overflow soil surface, mm
    % ----------

    % Constrains of temperature
    [s_tem] = temp_stress(Topt,Ta);

    % Case 0: groundwater overflow
    if zgw <= 0
        [wa, zgw, Tr, Es, uex] = swb_case0(wa, IWS, pEc, pEs, s_tem, soilpar, pftpar, wet, zm, zgw);
    
    % Case 1: groundwater table above layer 1
    elseif zgw>0 && zgw<= zm(1)
        [wa, zgw, Tr, Es, uex] = swb_case1(wa, IWS, pEc, pEs, s_tem, soilpar, pftpar, wet, zm, zgw);

    % Case 2: groundwater table between layer 1 and layer 2
    elseif zgw > zm(1) && zgw <= zm(1) + zm(2)
        [wa, zgw, Tr, Es, uex] = swb_case2(wa, IWS, pEc, pEs, s_tem, soilpar, pftpar, wet, zm, zgw);

    % Case 3: groundwater table between layer 2 and layer 3
    elseif zgw > zm(1) + zm(2) && zgw < zm(1) + zm(2) + zm(3)
        [wa, zgw, Tr, Es, uex] = swb_case3(wa, IWS, pEc, pEs, s_tem, soilpar, pftpar, wet, zm, zgw);

    % Case 4: groundwater table below layer 3
    else
        [wa, zgw, Tr, Es, uex] = swb_case4(wa, IWS, pEc, pEs, s_tem, soilpar, pftpar, wet, zm, zgw);

    end

end
