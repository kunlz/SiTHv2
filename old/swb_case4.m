function [wa, zgw, Tr, Es, uex] = swb_case4(wa, IWS, pEc, pEs, s_tem, soilpar, pftpar, wet, zm, zgw)
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
    theta_fc = soilpar(5);

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
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, pftpar, wet, zm);

    % Calculate the moisture constrains for plant and soil in unsaturated zone
    [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar);
    [f_sm2, ~] = swc_stress(wa2, soilpar);
    [f_sm3, ~] = swc_stress(wa3, soilpar);

    % actual transpiration
    Tr1 = f_sm1 .* s_tem .* Tr_p1;
    Tr2 = f_sm2 .* s_tem .* Tr_p2;
    Tr3 = f_sm3 .* s_tem .* Tr_p3;
    Tr = Tr1 + Tr2 + Tr3;

    % actual soil evaporation
    Es = f_sm_s1 .* pEs; % Only consider about the first layer

    % -------------------------------------- %
    % soil water drainage (unsaturated zone) %
    % -------------------------------------- %
    % layer #1
    % update the soil moisture after ET, layer #1
    wa1 = (wa1 * d1 - Es - Tr1) / d1;

    % drainage from unsaturated zone
    f1 = min(ks, d1 * (wa1 - theta_fc)); 
    f1 = max(f1, 0);

    % update the soil moisture after drainage, layer #1
    wa1 = (wa1 * d1 - f1) / d1;

    % layer #2
    % update the soil moisture after ET, layer #2
    wa2 = (wa2 * d2 - Tr2) / d2;

    % drainage from unsaturated zone
    f2 = min(ks, d2 * (wa2 - theta_fc)); 
    f2 = max(f2, 0);

    % update the soil moisture after drainage, layer #2
    wa2 = (wa2 * d2 + f1 - f2) / d2;

    if wa2 > theta_sat
        wa2 = theta_sat;
        ff2 = (wa2 - theta_sat) * d2;
    else
        ff2 = 0;
    end

    % layer #3
    % update the soil moisture after ET, layer #3, unsat-zone
    wa3 = (wa3 * d3 - Tr3) / d3;

    % drainage from unsaturated zone
    f3 = min(ks, d3 * (wa3 - theta_fc)); 
    f3 = max(f3, 0);

    % update the soil moisture after drainage, layer #3
    wa3 = (wa3 * d3 + f2 + ff2 - f3) / d3;

    if wa3 > theta_sat
        wa3 = theta_sat;
        ff3 = (wa3 - theta_sat) * d3;
    else
        ff3 = 0;
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

        wa3 = (wa3_unsat * (zgw - zm(1) - zm(2)) + theta_sat * (zm(1) + zm(2) + zm(3) - zgw)) / zm(3);

    elseif zgw > zm(1) && zgw < zm(1) + zm(2)

        wa2 = (wa2 * (zgw - zm(1)) + theta_sat * (zm(1) + zm(2) - zgw)) / zm(2);
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
