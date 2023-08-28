% ======================================================================= %
% ----------------------- %
%    Soil Water Balance   %     Kun Zhang, 2021/04/17, Beijing
% ----------------------- %
% ======================================================================= %
function [wa, zgw, Tr, Es, uex] = sw_balance2(IWS, pEc, pEs, Ta, Topt, s_VOD, ...
    wa, soilpar, pftpar, wet, zm, zgw)
    % ---------- function input -------
    % IWS     -- total water enter into soil surface, mm
    % pEc     -- potential ET allocate to plant, mm
    % pEs     -- potential ET allocate to soil surface, mm
    % Ta      -- air temperature, C
    % Topt    -- optimal growth temperature for plant, C
    % s_VOD   -- constrains of VOD,[0,1]
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
    [s_tem] = temp_stress(Topt, Ta);

    % Constrains of VPD
    % [s_vod] = vapor_stress(VPD, k2);
    % s_vod = 1;

    % Case 0: groundwater overflow
    if zgw <= 0
        [wa, zgw, Tr, Es, uex] = swb_case0(wa, IWS, pEc, pEs, s_tem, s_VOD, ...
            soilpar, pftpar, wet, zm, zgw);

        % Case 1: groundwater table in layer 1
    elseif zgw > 0 && zgw <= zm(1)
        [wa, zgw, Tr, Es, uex] = swb_case1(wa, IWS, pEc, pEs, s_tem, s_VOD, ...
            soilpar, pftpar, wet, zm, zgw);

        % Case 2: groundwater table in layer 2
    elseif zgw > zm(1) && zgw <= zm(1) + zm(2)
        [wa, zgw, Tr, Es, uex] = swb_case2(wa, IWS, pEc, pEs, s_tem, s_VOD, ...
            soilpar, pftpar, wet, zm, zgw);

        % Case 3: groundwater table in layer 3
    elseif zgw > zm(1) + zm(2) && zgw < zm(1) + zm(2) + zm(3)
        [wa, zgw, Tr, Es, uex] = swb_case3(wa, IWS, pEc, pEs, s_tem, s_VOD, ...
            soilpar, pftpar, wet, zm, zgw);

        % Case 4: groundwater table below layer 3
    else
        [wa, zgw, Tr, Es, uex] = swb_case4(wa, IWS, pEc, pEs, s_tem, s_VOD, ...
            soilpar, pftpar, wet, zm, zgw);

    end

end

% ======================================================================= %
% ------------------------------ %
% Case 0 -- groundwater overflow %
% ------------------------------ %
function [wa, zgw, Tr, Es, uex] = swb_case0(wa, IWS, pEc, pEs, s_tem, s_vod, ...
    soilpar, pftpar, wet, zm, zgw)
    % function input:
    % ----------
    % wa      -- soil water content, 3 layers
    % IWS     -- total water enter into soil surface, mm
    % pEc     -- potential ET allocate to plant, mm
    % pEs     -- potential ET allocate to soil surface, mm
    % soilpar -- soil-related parameters
    % pftpar  -- plant-related parameters
    % wet     -- wetness indice
    % zm      -- soil layer depth, 3 layers
    % zgw     -- groundwater table depth, mm
    % ----------

    % all saturated in layer #1 #2 #3

    % old soil water content in layer #1
    wa1 = wa(1);
    % old soil water content in layer #2
    wa2 = wa(2);
    % old soil water content in layer #3
    wa3 = wa(3);

    % hydraulic conductivity for specific soil type
    % ks = soilpar(1);

    % saturated swc for specific soil type
    theta_sat = soilpar(3);

    % field water capacity for specific soil type
    theta_fc = soilpar(5);

    % wilting point for specific soil type
    % wwp = soilpar(7);

    % ====== water supplement ====== %

    % layer #1 - saturated
    % full filled with groundwater

    % layer #2 - saturated
    % full filled with groundwater

    % layer #3 - saturated
    % full filled with groundwater

    % ====== water consumption ====== %

    % ------------------ %
    % Evapotranspiration %
    % ------------------ %

    % distributed the potential Tr to different layers
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, ...
    pftpar, wet, zm);

    % actual transpiration
    Tr1 = s_vod * s_tem * Tr_p1;
    Tr2 = s_vod * s_tem * Tr_p2;
    Tr3 = s_vod * s_tem * Tr_p3;
    Tr = Tr1 + Tr2 + Tr3;
    % actual soil evaporation
    Es = pEs; % only consider about the first layer

    % layer #1   full filled with groundwater
    % layer #2   full filled with groundwater
    % layer #3   full filled with groundwater

    % --------------------------- %
    % The groundwater table depth %
    % --------------------------- %

    % total water recharge to groundwater from unsaturated zone
    F1 = IWS;

    % total transpiration from groundwater
    Tr_g = Tr;

    % R_sb groundwater discaharge
    R_sb_max = 39; % mm day-1
    f = 1.25e-3; % mm-1
    R_sb = R_sb_max * exp(-f * zgw);

    % variation of water storaged in the saturated zone
    delta_w = F1 - Tr_g - R_sb; % should below zero

    % changes of groundwater table depth
    delta_zgw = delta_w / (theta_sat - theta_fc);
    zgw = zgw - delta_zgw;
    uex = 0; % excess water to soil surface, mm

    % update soil moisture and groundwater table depth
    if zgw > zm(1) + zm(2) + zm(3)

        wa1 = theta_fc;
        wa2 = theta_fc;
        wa3 = theta_fc;

    elseif zgw > zm(1) + zm(2) && zgw < zm(1) + zm(2) + zm(3)

        wa1 = theta_fc;
        wa2 = theta_fc;
        wa3 = (theta_fc * (zgw - zm(1) - zm(2)) + theta_sat * ...
            (zm(1) + zm(2) + zm(3) - zgw)) / zm(3);

    elseif zgw > zm(1) && zgw < zm(1) + zm(2)

        wa1 = theta_fc;
        wa2 = (theta_fc * (zgw - zm(1)) + theta_sat * (zm(1) + ...
            zm(2) - zgw)) / zm(2);
        wa3 = theta_sat;

    elseif zgw > 0 && zgw < zm(1)

        wa1 = (theta_fc * zgw + theta_sat * (zm(1) - zgw)) / zm(1);
        wa2 = theta_sat;
        wa3 = theta_sat;

    elseif zgw <= 0

        wa1 = theta_sat;
        wa2 = theta_sat;
        wa3 = theta_sat;

        uex = -zgw * theta_fc; % excess water to soil surface, mm

    end

    % updated soil water content
    wa = [wa1, wa2, wa3];
    zgw = max(0, zgw);

end

% ======================================================================= %
% --------------------------------------- %
% Case 1 -- groundwater table in layer 1  %
% --------------------------------------- %
function [wa, zgw, Tr, Es, uex] = swb_case1(wa, IWS, pEc, pEs, s_tem, s_vod, ...
    soilpar, pftpar, wet, zm, zgw)
    % function input:
    % ----------
    % wa      -- soil water content, 3 layers
    % IWS     -- total water enter into soil surface, mm
    % pEc     -- potential ET allocate to plant, mm
    % pEs     -- potential ET allocate to soil surface, mm
    % soilpar -- soil-related parameters
    % pftpar  -- plant-related parameters
    % wet     -- wetness indice
    % zm      -- soil layer depth, 3 layers
    % zgw     -- groundwater table depth, mm
    % ----------

    % unsaturated depth in layer #1
    d1 = zgw;

    % old soil water content in layer #1
    wa1 = wa(1);
    % old soil water content in layer #2
    wa2 = wa(2);
    % old soil water content in layer #3
    wa3 = wa(3);

    % hydraulic conductivity for specific soil type
    ks = soilpar(1);

    % saturated swc for specific soil type
    theta_sat = soilpar(3);

    % field water capacity for specific soil type
    theta_fc = soilpar(5);

    % wilting point for specific soil type
    wwp = soilpar(7);

    % ====== water supplement ====== %

    % layer #1
    % existed water column in the unsaturated zone #1
    wa1_unsat = (wa1 * zm(1) - theta_sat * (zm(1) - d1)) / d1;
    wc_s1 = d1 * wa1_unsat;

    % maximum water column in d1
    wc_m1 = d1 * theta_sat;

    if wc_s1 + IWS >= wc_m1

        % current soil water content
        wa1_unsat = theta_sat;

        % exceeded water
        vw1 = wc_s1 + IWS - wc_m1;
    else

        % soil water content in unsaturated zone
        wa1_unsat = wa1_unsat + IWS / d1;
        % calculate the adjusted swc#1 with considering the groundwater depth
        wa1 = (wa1_unsat * d1 + theta_sat * (zm(1) - d1)) / zm(1);
        % no exceeded water
        vw1 = 0;
    end

    % layer #2 - saturated
    % full filled with groundwater

    % layer #3 - saturated
    % full filled with groundwater

    % ====== water consumption ====== %

    % ------------------ %
    % Evapotranspiration %
    % ------------------ %

    % distributed the potential Tr to different layers
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, ...
    pftpar, wet, zm);

    % divide Tr_p1 into unsaturated zone and saturated zone
    Tr_p1_u = Tr_p1 * (d1 * wa1_unsat) / (d1 * wa1_unsat + (zm(1) - d1) ...
    * theta_sat);
    Tr_p1_g = Tr_p1 * ((zm(1) - d1) * theta_sat) / (d1 * wa1_unsat + ...
        (zm(1) - d1) * theta_sat);

    % calculate the moisture constrains for plant and soil in unsaturated zone
    [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar, pEc, pftpar);

    % actual transpiration
    Tr1_u = f_sm1 * s_vod * s_tem * Tr_p1_u;
    Tr1_u = max(Tr1_u, 0); % reject negtive value
    Tr1_u = min(Tr1_u, d1 * (wa1_unsat - wwp)); % less than maximum avaliable water
    Tr1_g = s_vod * s_tem * Tr_p1_g;
    Tr1 = Tr1_u + Tr1_g;
    Tr2 = s_vod * s_tem * Tr_p2;
    Tr3 = s_vod * s_tem * Tr_p3;

    Tr = Tr1 + Tr2 + Tr3;

    % actual soil evaporation
    Es = f_sm_s1 .* pEs; % only consider about the first layer
    Es_u = Es * (d1 * wa1_unsat) / (d1 * wa1_unsat + (zm(1) - d1) ...
        * theta_sat);
    Es_u = min(Es_u, d1 * wa1_unsat);
    Es_u = max(Es_u, 0);

    % -------------------------------------- %
    % soil water drainage (unsaturated zone) %
    % -------------------------------------- %

    % ---------------------------------------------------------------- layer #1
    % Ducharne % Polcher, 1998
    dd = 1.5;
    Dmin = 0.048; % mm day-1
    Dmax = 4.8; % mm day-1

    thx1 = wa1_unsat / theta_sat;

    if thx1 < 0.75
        Perc1 = Dmin * thx1;
    else
        Perc1 = Dmin * thx1 + (Dmax - Dmin) * thx1^dd;
    end

    % drainage from unsaturated zone, #1
    f1 = min(ks, Perc1);

    % update the soil moisture after ET & drainage, layer #1
    wa1_unsat = (wa1_unsat * d1 - f1 - Es_u - Tr1_u) / d1;
    wa1_unsat = max(wa1_unsat, 0);

    % layer #2   full filled with groundwater
    % layer #3   full filled with groundwater

    % --------------------------- %
    % The groundwater table depth %
    % --------------------------- %

    % total water recharge to groundwater
    F1 = f1 + vw1;

    % total transpiration from groundwater
    Tr_g = Tr1_g + Tr2 + Tr3;

    % R_sb groundwater discaharge
    R_sb_max = 39; % mm day-1
    f = 1.25e-3; % mm-1
    R_sb = R_sb_max * exp(-f * zgw);

    % variation of water storaged in the saturated zone
    delta_w = F1 - Tr_g - R_sb;

    % changes of groundwater table depth
%     if delta_w < 0 % decline of the water table
% 
%         delta_zgw = delta_w / theta_fc;
% 
%     else % increase of the water table

      delta_zgw = delta_w / (theta_sat - wa1_unsat);

%     end

    zgw = zgw - delta_zgw;
    uex = 0; % excess water to soil surface, mm

    % update soil moisture and groundwater table depth
    if zgw > zm(1) + zm(2) + zm(3)

        wa1 = (wa1_unsat * d1 + theta_fc * (zm(1) - d1)) / zm(1);
        wa2 = theta_fc;
        wa3 = theta_fc;

    elseif zgw > zm(1) + zm(2) && zgw < zm(1) + zm(2) + zm(3)

        wa1 = (wa1_unsat * d1 + theta_fc * (zm(1) - d1)) / zm(1);
        wa2 = theta_fc;
        wa3 = (theta_fc * (zgw - zm(1) - zm(2)) + theta_sat * (zm(1) + ...
            zm(2) + zm(3) - zgw)) / zm(3);

    elseif zgw > zm(1) && zgw < zm(1) + zm(2)

        wa1 = (wa1_unsat * d1 + theta_fc * (zm(1) - d1)) / zm(1);
        wa2 = (theta_fc * (zgw - zm(1)) + theta_sat * ...
            (zm(1) + zm(2) - zgw)) / zm(2);
        wa3 = theta_sat;

    elseif zgw > 0 && zgw < zm(1)

        wa1 = (wa1_unsat * zgw + theta_sat * (zm(1) - zgw)) / zm(1);
        wa2 = theta_sat;
        wa3 = theta_sat;

    elseif zgw <= 0

        wa1 = theta_sat;
        wa2 = theta_sat;
        wa3 = theta_sat;

        uex = -zgw * theta_fc; % excess water to soil surface, mm

    end

    % updated soil water content
    wa = [wa1, wa2, wa3];
    zgw = max(0, zgw);

end

% ======================================================================= %
% --------------------------------------- %
% Case 2 -- groundwater table in layer 2  %
% --------------------------------------- %
function [wa, zgw, Tr, Es, uex] = swb_case2(wa, IWS, pEc, pEs, s_tem, s_vod, ...
    soilpar, pftpar, wet, zm, zgw)
    % function input:
    % ----------
    % wa      -- soil water content, 3 layers
    % IWS     -- total water enter into soil surface, mm
    % pEc     -- potential ET allocate to plant, mm
    % pEs     -- potential ET allocate to soil surface, mm
    % soilpar -- soil-related parameters
    % pftpar  -- plant-related parameters
    % wet     -- wetness indice
    % zm      -- soil layer depth, 3 layers
    % zgw     -- groundwater table depth, mm
    % ----------

    % unsaturated depth in layer #1
    d1 = zm(1);
    % unsaturated depth in layer #2
    d2 = zgw - d1;

    % old soil water content in layer #1
    wa1 = wa(1);
    % old soil water content in layer #2
    wa2 = wa(2);
    % old soil water content in layer #3
    wa3 = wa(3);

    % hydraulic conductivity for specific soil type
    ks = soilpar(1);

    % saturated swc for specific soil type
    theta_sat = soilpar(3);

    % field water capacity for specific soil type
    theta_fc = soilpar(5);

    % wilting point for specific soil type
    wwp = soilpar(7);

    % ====== water supplement ====== %

    % layer #1
    % existed water column in the unsaturated zone #1
    wa1_unsat = wa1;
    wc_s1 = d1 * wa1_unsat;

    % maximum water column in d1
    wc_m1 = d1 * theta_sat;

    if wc_s1 + IWS >= wc_m1

        % current soil water content
        wa1 = theta_sat;
        % exceeded water
        vw1 = wc_s1 + IWS - wc_m1;
    else

        % soil water content in unsaturated zone
        wa1 = wa1_unsat + IWS / d1;
        % no exceeded water
        vw1 = 0;
    end
    
    % vw1 = vw1 + IWU;
    % layer #2
    % existed water column in the unsaturated zone #2
    wa2_unsat = (wa2 * zm(2) - theta_sat * (zm(2) - d2)) / d2;
    wc_s2 = d2 * wa2_unsat;

    % maximum water column in d2
    wc_m2 = d2 * theta_sat;

    if wc_s2 + vw1 >= wc_m2

        % current soil water content
        wa2_unsat = theta_sat;
        % exceeded water
        vw2 = wc_s2 + vw1 - wc_m2;
    else

        % soil water content in unsaturated zone
        wa2_unsat = wa2_unsat + vw1 / d2;
        % calculate the adjusted swc#2 with considering the groundwater depth
        wa2 = (wa2_unsat * d2 + theta_sat * (zm(2) - d2)) / zm(2);
        % no exceeded water
        vw2 = 0;
    end

    % layer #3
    % full filled with groundwater

    % ====== water consumption ====== %

    % ------------------ %
    % Evapotranspiration %
    % ------------------ %

    % distributed the potential Tr to different layers
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, ...
    pftpar, wet, zm);

    % divide Tr_p2 into unsaturated zone and saturated zone
    Tr_p2_u = Tr_p2 * (d2 * wa2_unsat) / (d2 * wa2_unsat + ...
    (zm(2) - d2) * theta_sat);
    Tr_p2_g = Tr_p2 * ((zm(2) - d2) * theta_sat) / (d2 * wa2_unsat + ...
        (zm(2) - d2) * theta_sat);

    % calculate the moisture constrains for plant and soil in unsaturated zone
    [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar, pEc, pftpar);
    [f_sm2, ~] = swc_stress(wa2_unsat, soilpar, pEc, pftpar);

    % actual transpiration
    Tr1 = f_sm1 * s_vod * s_tem * Tr_p1;
    Tr2_u = f_sm2 * s_vod * s_tem * Tr_p2_u;
    Tr2_g = s_vod * s_tem * Tr_p2_g;
    Tr2 = Tr2_u + Tr2_g;
    Tr3 = s_vod * s_tem * Tr_p3;
    Tr = Tr1 + Tr2 + Tr3;

    % actual soil evaporation
    Es = f_sm_s1 .* pEs; % only considering the first layer

    % -------------------------------------- %
    % soil water drainage (unsaturated zone) %
    % -------------------------------------- %

    % ---------------------------------------------------------------- layer #1
    % check Es & Tr, layer #1   unsat-zone
    Es = max(Es, 0);
    Tr1 = max(Tr1, 0);

    if wa1 > 0 && Es + Tr1 > d1 * wa1  % wilting point (0) for first layer
        Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es);
        Es = d1 * wa1 - Tr1;
    end

    % Ducharne % Polcher, 1998
    dd = 1.5;
    Dmin = 0.048; % mm day-1
    Dmax = 4.8; % mm day-1

    thx1 = wa1 / theta_sat;

    if thx1 < 0.75
        Perc1 = Dmin * thx1;
    else
        Perc1 = Dmin * thx1 + (Dmax - Dmin) * thx1^dd;
    end

    % drainage from unsaturated zone, #1
    f1 = min(ks, Perc1);

    % update the soil moisture after drainage, layer #1
    wa1 = (wa1 * d1 - f1 - Es - Tr1) / d1;
    wa1 = max(wa1, 0);

    % ---------------------------------------------------------------- layer #2
    % check Tr, layer #2   unsat-zone
    Tr2_u = max(Tr2_u, 0); % reject negtive value
    Tr2_u = min(Tr2_u, d2 * (wa2_unsat - wwp)); % less than maximum avaliable water

    % gravity drainage
    Dmin = 0.012; % 0.0005*24, mm day-1
    Dmax = 1.2; % 0.05*24,   mm day-1
    thx2 = wa2_unsat / theta_sat;

    if thx2 < 0.75
        Perc2 = Dmin * thx2;
    else
        Perc2 = Dmin * thx2 + (Dmax - Dmin) * thx2^dd;
    end

    % drainage from unsaturated zone, #2
    f2 = min(ks, Perc2);

    % update the soil moisture after ET & drainage, layer #2 un
    wa2_unsat = (wa2_unsat * d2 + f1 - f2 - Tr2_u) / d2;

    wa2_unsat = max(wa2_unsat, 0); % > 0
    wa2_unsat = min(wa2_unsat, 1); % < 1

    if wa2_unsat > theta_sat
        ff2 = (wa2_unsat - theta_sat) * d2; % extra water from upper layer
        ff2 = max(ff2, 0);
        wa2_unsat = theta_sat;
    else
        ff2 = 0;
    end

    % layer #3
    % full filled with groundwater

    % --------------------------- %
    % The groundwater table depth %
    % --------------------------- %

    % total water recharge to groundwater
    F1 = f2 + ff2 + vw2;

    % total transpiration from groundwater
    Tr_g = Tr2_g + Tr3;

    % R_sb groundwater discaharge
    R_sb_max = 39; % mm day-1
    f = 1.25e-3; % mm-1
    R_sb = R_sb_max * exp(-f * zgw);

    % variation of water storaged in the saturated zone
    delta_w = F1 - Tr_g - R_sb;

    % changes of groundwater table depth
%     if delta_w < 0 % decline of the water table
% 
%         delta_zgw = delta_w / (theta_sat - theta_fc);
% 
%     else % increase of the water table

    delta_zgw = delta_w / (theta_sat - (wa1 + wa2_unsat) / 2);

%     end

    zgw = zgw - delta_zgw;
    uex = 0; % excess water to soil surface, mm

    % update soil moisture and groundwater table depth
    if zgw > zm(1) + zm(2) + zm(3)

        wa2 = (wa2_unsat * d2 + theta_fc * (zm(2) - d2)) / zm(2);
        wa3 = theta_fc;

    elseif zgw > zm(1) + zm(2) && zgw < zm(1) + zm(2) + zm(3)

        wa2 = (wa2_unsat * d2 + theta_fc * (zm(2) - d2)) / zm(2);
        wa3 = (theta_fc * (zgw - zm(1) - zm(2)) + theta_sat * ...
            (zm(1) + zm(2) + zm(3) - zgw)) / zm(3);

    elseif zgw > zm(1) && zgw < zm(1) + zm(2)

        wa2 = (wa2_unsat * (zgw - zm(1)) + theta_sat * ...
            (zm(1) + zm(2) - zgw)) / zm(2);
        wa3 = theta_sat;

    elseif zgw > 0 && zgw < zm(1)

        wa1 = (wa1 * zgw + theta_sat * (zm(1) - zgw)) / zm(1);
        wa2 = theta_sat;
        wa3 = theta_sat;

    elseif zgw <= 0

        wa1 = theta_sat;
        wa2 = theta_sat;
        wa3 = theta_sat;

        uex = -zgw * theta_fc; % excess water to soil surface, mm

    end

    % updated soil water content
    wa = [wa1, wa2, wa3];
    zgw = max(0, zgw);

end

% ======================================================================= %
% --------------------------------------- %
% Case 3 -- groundwater table in layer 3  %
% --------------------------------------- %
function [wa, zgw, Tr, Es, uex] = swb_case3(wa, IWS, pEc, pEs, s_tem, s_vod, ...
    soilpar, pftpar, wet, zm, zgw)
    % function input:
    % ----------
    % wa      -- soil water content, 3 layers
    % IWS     -- total water enter into soil surface, mm
    % pEc     -- potential ET allocate to plant, mm
    % pEs     -- potential ET allocate to soil surface, mm
    % soilpar -- soil-related parameters
    % pftpar  -- plant-related parameters
    % wet     -- wetness indice
    % zm      -- soil layer depth, 3 layers
    % zgw     -- groundwater table depth, mm
    % ----------

    % unsaturated depth in layer #1
    d1 = zm(1);
    % unsaturated depth in layer #2
    d2 = zm(2);
    % unsaturated depth in layer #3
    d3 = zgw - zm(1) - zm(2);

    % old soil water content in layer #1
    wa1 = wa(1);
    % old soil water content in layer #2
    wa2 = wa(2);
    % old soil water content in layer #3
    wa3 = wa(3);

    % hydraulic conductivity for specific soil type
    ks = soilpar(1);

    % saturated swc for specific soil type
    theta_sat = soilpar(3);

    % field water capacity for specific soil type
    theta_fc = soilpar(5);

    % wilting point for specific soil type
    wwp = soilpar(7);

    % ====== water supplement ====== %

    % layer #1
    % existed water column in the unsaturated zone #1
    wa1_unsat = wa1;
    wc_s1 = d1 * wa1_unsat;

    % maximum water column in d1
    wc_m1 = d1 * theta_sat;

    if wc_s1 + IWS >= wc_m1

        % current soil water content
        wa1 = theta_sat;
        % exceeded water
        vw1 = wc_s1 + IWS - wc_m1;
    else

        % soil water content in unsaturated zone
        wa1 = wa1_unsat + IWS / d1;
        vw1 = 0; % no exceeded water
    end
    
    % vw1 = vw1+IWU;
    % layer #2
    % existed water column in the unsaturated zone #2
    wa2_unsat = wa2;
    wc_s2 = d2 * wa2_unsat;

    % maximum water column in d2
    wc_m2 = d2 * theta_sat;

    if wc_s2 + vw1 >= wc_m2

        % current soil water content
        wa2 = theta_sat;
        % exceeded water
        vw2 = wc_s2 + vw1 - wc_m2;
    else
        % soil water content in unsaturated zone
        wa2 = wa2_unsat + vw1 / d2;
        vw2 = 0; % no exceeded water
    end

    % layer #3
    % existed water column in the unsaturated zone #3
    wa3_unsat = (wa3 * zm(3) - theta_sat * (zm(3) - d3)) / d3;
    wc_s3 = d3 * wa3_unsat;

    % maximum water column in d3
    wc_m3 = d3 * theta_sat;

    if wc_s3 + vw2 >= wc_m3
        wa3 = theta_sat;
        vw3 = wc_s3 + vw2 - wc_m3;
    else
        % soil water content in unsaturated zone
        wa3_unsat = wa3_unsat + vw2 / d3;
        % calculate the adjusted swc#3 with considering the groundwater depth
        wa3 = (wa3_unsat * d3 + theta_sat * (zm(3) - d3)) / zm(3);
        % no exceeded water
        vw3 = 0;
    end

    % ====== water consumption ====== %

    % ------------------ %
    % Evapotranspiration %
    % ------------------ %

    % distributed the potential Tr to different layers
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, ...
    pftpar, wet, zm);

    % divide Tr_p3 into unsaturated and saturated zone at the layer #3
    Tr_p3_u = Tr_p3 * (d3 * wa3_unsat) / (d3 * wa3_unsat + ...
    (zm(3) - d3) * theta_sat);
    Tr_p3_g = Tr_p3 * ((zm(3) - d3) * theta_sat) / (d3 * wa3_unsat + ...
        (zm(3) - d3) * theta_sat);

    % calculate the moisture constrains for plant and soil in unsaturated zone
    [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar, pEc, pftpar);
    [f_sm2, ~] = swc_stress(wa2, soilpar, pEc, pftpar);
    [f_sm3, ~] = swc_stress(wa3_unsat, soilpar, pEc, pftpar);

    % actual transpiration
    Tr1 = f_sm1 * s_vod * s_tem * Tr_p1;
    Tr2 = f_sm2 * s_vod * s_tem * Tr_p2;
    Tr3_u = f_sm3 * s_vod * s_tem * Tr_p3_u;
    Tr3_g = s_vod * s_tem * Tr_p3_g;
    Tr3 = Tr3_u + Tr3_g;
    Tr = Tr1 + Tr2 + Tr3;

    % actual soil evaporation
    Es = f_sm_s1 .* pEs; % only considering the first layer

    % -------------------------------------- %
    % soil water drainage (unsaturated zone) %
    % -------------------------------------- %
    % ---------------------------------------------------------------- layer #1
    Es = max(Es, 0);
    Tr1 = max(Tr1, 0);

    if wa1 > 0 && Es + Tr1 > d1 * wa1 % wilting point -- 0, first layer
        Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es);
        Es = d1 * wa1 - Tr1;
    end

    % Ducharne & Polcher, 1998
    dd = 1.5;
    Dmin = 0.048; % mm day-1
    Dmax = 4.8; % mm day-1

    thx1 = wa1 / theta_sat;

    if thx1 < 0.75
        Perc1 = Dmin * thx1;
    else
        Perc1 = Dmin * thx1 + (Dmax - Dmin) * thx1^dd;
    end

    % drainage from unsaturated zone, #1
    f1 = min(ks, Perc1);

    % update the soil moisture after ET & drainage, layer #1
    wa1 = (wa1 * d1 - f1 - Es - Tr1) / d1;
    wa1 = max(wa1, 0);

    % ---------------------------------------------------------------- layer #2
    Tr2 = max(Tr2, 0); % reject negtive value
    Tr2 = min(Tr2, d2 * (wa2 - wwp)); % less than maximum avaliable water

    % gravity drainage
    Dmin = 0.012; % 0.0005*24, mm day-1
    Dmax = 1.2; % 0.05*24,   mm day-1
    thx2 = wa2 / theta_sat;

    if thx2 < 0.75
        Perc2 = Dmin * thx2;
    else
        Perc2 = Dmin * thx2 + (Dmax - Dmin) * thx2^dd;
    end

    % drainage from unsaturated zone, #2
    f2 = min(ks, Perc2);

    % update the soil moisture after ET & drainage, layer #2
    wa2 = (wa2 * d2 + f1 - f2 - Tr2) / d2;

    wa2 = max(wa2, 0); % > wwp--0
    wa2 = min(wa2, 1); % < 1

    if wa2 > theta_sat
        ff2 = (wa2 - theta_sat) * d2; % extra water from upper layer
        ff2 = max(ff2, 0);
        wa2 = theta_sat;
    else
        ff2 = 0;
    end

    % ---------------------------------------------------------------- layer #3
    % check Tr, layer #3   unsat-zone
    Tr3_u = max(Tr3_u, 0); % reject negtive value
    Tr3_u = min(Tr3_u, d3 * (wa3_unsat - wwp)); % less than maximum avaliable water

    % gravity drainage
    Dmin = 0.012; % 0.0005*24, mm day-1
    Dmax = 1.2; % 0.05*24,   mm day-1
    thx3 = wa3_unsat / theta_sat;

    if thx3 < 0.75
        Perc3 = Dmin * thx3;
    else
        Perc3 = Dmin * thx3 + (Dmax - Dmin) * thx3^dd;
    end

    % drainage from unsaturated zone, #3
    f3 = min(ks, Perc3);

    % update the soil moisture after ET & drainage, layer #3
    wa3_unsat = (wa3_unsat * d3 + f2 + ff2 - f3 - Tr3_u) / d3;

    wa3_unsat = max(wa3_unsat, 0); % > wwp 0
    wa3_unsat = min(wa3_unsat, 1); % < theta_s 1
    
    if wa3_unsat > theta_sat
        ff3 = (wa3_unsat - theta_sat) * d3; % extra water from upper layer
        ff3 = max(ff3, 0);
        wa3_unsat = theta_sat;
    else
        ff3 = 0;
    end

    % --------------------------- %
    % The groundwater table depth %
    % --------------------------- %

    % total water recharge to groundwater
    F1 = f3 + ff3 + vw3;

    % total transpiration from groundwater
    Tr_g = Tr3_g;

    % R_sb groundwater discaharge
    R_sb_max = 39; % mm day-1
    f = 1.25e-3; % mm-1
    R_sb = R_sb_max * exp(-f * zgw);

    % variation of water storaged in the saturated zone
    delta_w = F1 - Tr_g - R_sb;

    % changes of groundwater table depth
%     if delta_w < 0 % decline of the water table
% 
%         delta_zgw = delta_w / (theta_sat - theta_fc);
% 
%     else % increase of the water table

    delta_zgw = delta_w / (theta_sat - (wa1 + wa2 + wa3_unsat) / 3);

%     end

    zgw = zgw - delta_zgw;
    uex = 0; % excess water to soil surface, mm

    % update soil moisture and groundwater table depth
    if zgw > zm(1) + zm(2) + zm(3)

        wa3 = (wa3_unsat * d3 + theta_fc * (zm(3) - d3)) / zm(3);

    elseif zgw > zm(1) + zm(2) && zgw < zm(1) + zm(2) + zm(3)

        wa3 = (wa3_unsat * (zgw - zm(1) - zm(2)) + theta_sat * ...
            (zm(1) + zm(2) + zm(3) - zgw)) / zm(3);

    elseif zgw > zm(1) && zgw < zm(1) + zm(2)

        wa2 = (wa2 * (zgw - zm(1)) + theta_sat * (zm(1) + zm(2) - zgw)) ...
            / zm(2);
        wa3 = theta_sat;

    elseif zgw > 0 && zgw < zm(1)

        wa1 = (wa1 * zgw + theta_sat * (zm(1) - zgw)) / zm(1);
        wa2 = theta_sat;
        wa3 = theta_sat;

    elseif zgw <= 0

        wa1 = theta_sat;
        wa2 = theta_sat;
        wa3 = theta_sat;

        uex = -zgw * theta_sat; % excess water to soil surface, mm

    end

    % updated soil water content
    wa = [wa1, wa2, wa3];
    zgw = max(0, zgw);

end

% ======================================================================= %
% --------------------------------------- %
% Case 4 -- groundwater table in layer 4  %
% --------------------------------------- %
function [wa, zgw, Tr, Es, uex] = swb_case4(wa, IWS, pEc, pEs, s_tem, s_vod, ...
    soilpar, pftpar, wet, zm, zgw)
    % function input:
    % ----------
    % wa      -- soil water content, 3 layers
    % IWS     -- total water enter into soil surface, mm
    % pEc     -- potential ET allocate to plant, mm
    % pEs     -- potential ET allocate to soil surface, mm
    % soilpar -- soil-related parameters
    % pftpar  -- plant-related parameters
    % wet     -- wetness indice
    % zm      -- soil layer depth, 3 layers
    % zgw     -- groundwater table depth, mm
    % ----------

    % unsaturated depth in layer #1
    d1 = zm(1);
    % unsaturated depth in layer #2
    d2 = zm(2);
    % unsaturated depth in layer #3
    d3 = zm(3);

    % old soil water content in layer #1
    wa1 = wa(1);
    % old soil water content in layer #2
    wa2 = wa(2);
    % old soil water content in layer #3
    wa3 = wa(3);

    % hydraulic conductivity for specific soil type
    ks = soilpar(1);

    % saturated swc for specific soil type
    theta_sat = soilpar(3);

    % field water capacity for specific soil type
    % theta_fc = soilpar(5);

    % wilting point for specific soil type
    wwp = soilpar(7);

    % ====== water supplement ====== %

    % layer #1
    % existed water column in the unsaturated zone #1
    wa1_unsat = wa1;
    wc_s1 = d1 * wa1_unsat;

    % maximum water column in d1
    wc_m1 = d1 * theta_sat;

    if wc_s1 + IWS >= wc_m1 % saturated

        % current soil water content
        wa1 = theta_sat;
        % exceeded water
        vw1 = wc_s1 + IWS - wc_m1;

    else % un-saturated

        % soil water content in unsaturated zone
        wa1 = wa1_unsat + IWS / d1;
        % no exceeded water
        vw1 = 0;
    end

    % vw1 = vw1+IWU;
    % layer #2
    % existed water column in the unsaturated zone #2
    wa2_unsat = wa2;
    wc_s2 = d2 * wa2_unsat;

    % maximum water column in d2
    wc_m2 = d2 * theta_sat;

    if wc_s2 + vw1 >= wc_m2

        % current soil water content
        wa2 = theta_sat;
        % exceeded water to layer #3
        vw2 = wc_s2 + vw1 - wc_m2;
    else

        % soil water content in unsaturated zone
        wa2 = wa2_unsat + vw1 / d2;
        % no exceeded water
        vw2 = 0;
    end

    % layer #3
    % existed water column in the unsaturated zone #3
    wa3_unsat = wa3;
    wc_s3 = d3 * wa3_unsat;

    % maximum water column in d3
    wc_m3 = d3 * theta_sat;

    if wc_s3 + vw2 >= wc_m3

        % current soil water content
        wa3 = theta_sat;
        % exceeded water
        vw3 = wc_s3 + vw2 - wc_m3;
    else

        % soil water content in unsaturated zone
        wa3 = wa3_unsat + vw2 / d3;
        % no exceeded water
        vw3 = 0;
    end

    % ====== water consumption ====== %

    % ------------------ %
    % Evapotranspiration %
    % ------------------ %

    % distributed the potential T to different layers
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, ...
    pftpar, wet, zm);

    % Calculate the moisture constrains for plant and soil in unsaturated zone
    [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar, pEc, pftpar);
    [f_sm2, ~] = swc_stress(wa2, soilpar, pEc, pftpar);
    [f_sm3, ~] = swc_stress(wa3, soilpar, pEc, pftpar);

    % actual transpiration
    Tr1 = f_sm1 .* s_vod * s_tem .* Tr_p1;
    Tr2 = f_sm2 .* s_vod * s_tem .* Tr_p2;
    Tr3 = f_sm3 .* s_vod * s_tem .* Tr_p3;
    Tr = Tr1 + Tr2 + Tr3;

    % actual soil evaporation
    Es = f_sm_s1 .* pEs; % Only consider about the first layer

    % -------------------------------------- %
    % soil water drainage (unsaturated zone) %
    % -------------------------------------- %
    % ---------------------------------------------------------------- layer #1
    % update the soil moisture after ET, layer #1
    Es = max(Es, 0);
    Tr1 = max(Tr1, 0);

    if wa1 > 0 && Es + Tr1 > d1 * wa1 % wilting point
        Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es);
        Es = d1 * wa1 - Tr1;
    end

    % Ducharne % Polcher, 1998
    dd = 1.5;
    Dmin = 0.048; % mm day-1
    Dmax = 4.8; % mm day-1

    thx1 = wa1 / theta_sat;

    if thx1 < 0.75
        Perc1 = Dmin * thx1;
    else
        Perc1 = Dmin * thx1 + (Dmax - Dmin) * thx1^dd;
    end

    % drainage from unsaturated zone, #1
    f1 = min(ks, Perc1);

    % update the soil moisture after drainage, layer #1
    wa1 = (wa1 * d1 - f1 - Es - Tr1) / d1;
    wa1 = max(wa1, 0);
    wa1 = min(wa1, 1); % < theta_s

    % ---------------------------------------------------------------- layer #2
    % update the soil moisture after ET, layer #2
    Tr2 = max(Tr2, 0); % reject negtive value
    Tr2 = min(Tr2, d2 * (wa2 - wwp)); % less than maximum avaliable water

    % gravity drainage
    Dmin = 0.012; % 0.0005*24, mm day-1
    Dmax = 1.2; % 0.05*24,   mm day-1
    thx2 = wa2 / theta_sat;

    if thx2 < 0.75
        Perc2 = Dmin * thx2;
    else
        Perc2 = Dmin * thx2 + (Dmax - Dmin) * thx2^dd;
    end

    % drainage from unsaturated zone, #2
    f2 = min(ks, Perc2);

    % update the soil moisture after drainage, layer #2
    wa2 = (wa2 * d2 + f1 - f2 - Tr2) / d2;

    wa2 = max(wa2, 0); % > wwp 0
    wa2 = min(wa2, 1); % < theta_s 1

    if wa2 > theta_sat
        ff2 = (wa2 - theta_sat) * d2; % extra water from upper layer
        ff2 = max(ff2, 0);
        wa2 = theta_sat;
    else
        ff2 = 0;
    end

    % ---------------------------------------------------------------- layer #3
    % update the soil moisture after ET, layer #3, unsat-zone
    Tr3 = max(Tr3, 0); % reject negtive value
    Tr3 = min(Tr3, d3 * (wa3 - wwp)); % less than maximum avaliable water

    % gravity drainage
    Dmin = 0.012; % 0.0005*24, mm day-1
    Dmax = 1.2; % 0.05*24,   mm day-1
    thx3 = wa3 / theta_sat;

    if thx3 < 0.75
        Perc3 = Dmin * thx3;
    else
        Perc3 = Dmin * thx3 + (Dmax - Dmin) * thx3^dd;
    end

    % drainage from unsaturated zone, #3
    f3 = min(ks, Perc3);

    % update the soil moisture after drainage, layer #3
    wa3 = (wa3 * d3 + f2 + ff2 - f3 - Tr3) / d3;

    wa3 = max(wa3, 0); % > wwp
    wa3 = min(wa3, 1); % < theta_s
    if wa3 > theta_sat
        ff3 = (wa3 - theta_sat) * d3; % extra water from upper layer
        ff3 = max(ff3, 0);
        wa3_unsat = theta_sat;
    else
        ff3 = 0;
        wa3_unsat = wa3;
    end

    % --------------------------- %
    % The groundwater table depth %
    % --------------------------- %

    % total water recharge to groundwater
    F1 = f3 + ff3 + vw3;

    % total transpiration from groundwater
    Tr_g = 0;

    % R_sb groundwater discaharge
    R_sb_max = 39; % mm day-1
    f = 1.25e-3; % mm-1
    R_sb = R_sb_max * exp(-f * zgw);

    % variation of water storaged in the saturated zone
    delta_w = F1 - Tr_g - R_sb;

    % changes in groundwater table depth
    delta_zgw = delta_w / 0.2;
    zgw = zgw - delta_zgw;
    uex = 0; % excess water to soil surface, mm

    % update soil moisture and groundwater table depth
    if zgw > zm(1) + zm(2) && zgw < zm(1) + zm(2) + zm(3)

        wa3 = (wa3_unsat * (zgw - zm(1) - zm(2)) + theta_sat * ...
            (zm(1) + zm(2) + zm(3) - zgw)) / zm(3);

    elseif zgw > zm(1) && zgw < zm(1) + zm(2)

        wa2 = (wa2 * (zgw - zm(1)) + theta_sat * (zm(1) + zm(2) - zgw)) ...
            / zm(2);
        wa3 = theta_sat;

    elseif zgw > 0 && zgw < zm(1)

        wa1 = (wa1 * zgw + theta_sat * (zm(1) - zgw)) / zm(1);
        wa2 = theta_sat;
        wa3 = theta_sat;

    elseif zgw <= 0

        wa1 = theta_sat;
        wa2 = theta_sat;
        wa3 = theta_sat;

        uex = -zgw * theta_sat; % excess water to soil surface, mm

    end

    % updated soil water content
    wa = [wa1, wa2, wa3];
    zgw = max(0, zgw);

end

% ========================== E N D ====================================== %
