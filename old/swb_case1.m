function [wa, zgw, Tr, Es, uex] = swb_case1(wa, IWS, pEc, pEs, s_tem, soilpar, pftpar, wet, zm, zgw)
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
    [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, soilpar, pftpar, wet, zm);

    % divide Tr_p1 into unsaturated zone and saturated zone
    Tr_p1_u = Tr_p1 * (d1 * wa1_unsat) / (d1 * wa1_unsat + (zm(1) - d1) * theta_sat);
    Tr_p1_g = Tr_p1 * ((zm(1) - d1) * theta_sat) / (d1 * wa1_unsat + (zm(1) - d1) * theta_sat);

    % calculate the moisture constrains for plant and soil in unsaturated zone
    [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar);

    % actual transpiration
    Tr1_u = f_sm1 * s_tem * Tr_p1_u;
    Tr1_g = s_tem * Tr_p1_g;
    Tr1 = Tr1_u + Tr1_g;
    Tr2 = s_tem * Tr_p2;
    Tr3 = s_tem * Tr_p3;
    Tr = Tr1 + Tr2 + Tr3;

    % actual soil evaporation
    Es = f_sm_s1 .* pEs; % only consider about the first layer

    % -------------------------------------- %
    % soil water drainage (unsaturated zone) %
    % -------------------------------------- %

    % layer #1
    % update the soil moisture after ET, layer #1
    wa1_unsat = (wa1_unsat * d1 - Tr1 - Es) / d1;

    % drainage from unsaturated zone
    f1 = min(ks, d1 * (wa1_unsat - theta_fc));
    f1 = max(f1, 0);

    % update the soil moisture after drainage, layer #1
    wa1_unsat = (wa1_unsat * d1 - f1) / d1;

    if wa1_unsat > theta_sat
        wa1_unsat = theta_sat;
        ff1 = (wa1_unsat - theta_sat) * d1;
    else
        ff1 = 0;
    end

    % layer #2   full filled with groundwater
    % layer #3   full filled with groundwater

    % --------------------------- %
    % The groundwater table depth %
    % --------------------------- %

    % total water recharge to groundwater
    F1 = f1 + ff1 + vw1;

    % total transpiration from groundwater
    Tr_g = Tr1_g + Tr2 + Tr3;

    % R_sb groundwater discaharge
    R_sb_max = 39; % mm day-1
    f = 1.25e-3; % mm-1
    R_sb = R_sb_max * exp(-f * zgw);

    % variation of water storaged in the saturated zone
    delta_w = F1 - Tr_g - R_sb;

    % changes of groundwater table depth
    if delta_w < 0 % decline of the water table

        delta_zgw = delta_w / (theta_sat - theta_fc);

    else % increase of the water table

        delta_zgw = delta_w / (theta_sat - wa1_unsat);

    end
    zgw = zgw - delta_zgw; 
    uex = 0; % excess water to soil surface, mm

    % update soil moisture and groundwater table depth
    if zgw > zm(1) + zm(2) + zm(3)

        wa1 = (wa1_unsat * d1 + theta_fc * (zm(1)  - d1)) / zm(1);
        wa2 = theta_fc;
        wa3 = theta_fc;

    elseif zgw > zm(1) + zm(2) && zgw < zm(1) + zm(2) + zm(3)

        wa1 = (wa1_unsat * d1 + theta_fc * (zm(1)  - d1)) / zm(1);
        wa2 = theta_fc;
        wa3 = (theta_fc * (zgw - zm(1) - zm(2)) + theta_sat * (zm(1) + zm(2) + zm(3) - zgw)) / zm(3);

    elseif zgw > zm(1) && zgw < zm(1) + zm(2)

        wa1 = (wa1_unsat * d1 + theta_fc * (zm(1)  - d1)) / zm(1);
        wa2 = (theta_fc * (zgw - zm(1)) + theta_sat * (zm(1) + zm(2) - zgw)) / zm(2);
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
